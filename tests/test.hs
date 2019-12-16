{-# LANGUAGE RecordWildCards, ScopedTypeVariables, OverloadedStrings,
             ViewPatterns, ImplicitParams, OverloadedLists #-}
import Test.Tasty
import Test.Tasty.HUnit

import qualified Numeric.Sundials.ARKode.ODE as ARK
import qualified Numeric.Sundials.CVode.ODE  as CV
import Numeric.Sundials.Types
import Numeric.LinearAlgebra as L
import qualified Data.Vector.Storable as V
import Katip
import Foreign.C.Types
import System.IO
import Text.Printf (printf)
import GHC.Stack
import Control.Monad

----------------------------------------------------------------------
--                            Helpers
----------------------------------------------------------------------

data OdeProblem = OdeProblem
  { odeEvents :: [EventSpec]
  , odeMaxEvents :: !Int
  , odeRhs :: OdeRhs
  , odeJacobian :: Maybe (Double -> Vector Double -> Matrix Double)
  , odeInitCond :: V.Vector Double
  , odeSolTimes :: V.Vector Double
  , odeTolerances :: StepControl
  }

solveCV
  :: Katip m
  => ODEOpts CV.ODEMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCV opts OdeProblem{..} =
  CV.odeSolveWithEvents opts {stepControl = odeTolerances} odeEvents odeMaxEvents odeRhs odeJacobian odeInitCond odeSolTimes

main = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  log_env <- registerScribe "stderr" handleScribe defaultScribeSettings =<<
    initLogEnv "test" "devel"
  let ?log_env = log_env

  defaultMain $ testGroup "Tests" $
    [
      testGroup solver_name
      [ testGroup (show method) $
          let opts = defaultOpts method in
          [ withVsWithoutJacobian opts solver
          , eventTests opts solver
          , noErrorTests opts solver
          ]
      | method <- [CV.BDF, CV.ADAMS]
      ]
    | (solver, solver_name) <-
      [ (solveCV, "CVode") ]
    ] ++
    [ testGroup "Method comparison"
      [ testGroup (show method1 ++ " vs " ++ show method2) $ compareMethodsTests
          (defaultOpts method1, solver1)
          (defaultOpts method2, solver2)
      | (method1, solver1) <- availableMethods
      , (method2, solver2) <- availableMethods
      , method1 /= method2
      ]
    ]

availableMethods =
  [ (CV.BDF, solveCV)
  , (CV.ADAMS, solveCV)
  ]

defaultOpts :: method -> ODEOpts method
defaultOpts method = ODEOpts
  { maxNumSteps = 10000
  , minStep     = 1.0e-12
  , maxFail     = 10
  , odeMethod   = method
  , stepControl = defaultTolerances
  , initStep    = Nothing
  }

defaultTolerances = CV.XX' 1.0e-6 1.0e-10 1 1

checkDiscrepancy :: HasCallStack => Double -> Double -> Assertion
checkDiscrepancy eps diff = assertBool msg $ diff <= eps
  where
    msg = printf "Difference too large: %.2e > %.2e"
      diff eps

----------------------------------------------------------------------
--                             The tests
----------------------------------------------------------------------

noErrorTests opts solver = testGroup "Absence of error"
  [ testCase name $ do
      r <- runKatipT ?log_env $ solver opts prob
      case r of
        Right _ -> return ()
        Left e -> assertFailure (show e)
  | (name, prob) <- [ empty ]
  ]

withVsWithoutJacobian opts solver = testGroup "With vs without jacobian"
  [ testCase name $ do
      Right (solutionMatrix -> solJac)   <- runKatipT ?log_env $ solver opts prob
      Right (solutionMatrix -> solNoJac) <- runKatipT ?log_env $ solver opts prob { odeJacobian = Nothing }
      checkDiscrepancy 1e-3 $ norm_2 (solJac - solNoJac)
  | (name, prob) <- [ brusselator, robertson ]
  ]

compareMethodsTests (opts1, solver1) (opts2, solver2) =
  [ testCase name $ do
      Right (solutionMatrix -> sol1) <- runKatipT ?log_env $ solver1 opts1 prob
      Right (solutionMatrix -> sol2) <- runKatipT ?log_env $ solver2 opts2 prob
      let diff = maximum $ map abs $
                 zipWith (-) ((toLists $ tr sol1)!!0) ((toLists $ tr sol2)!!0)
      checkDiscrepancy 1e-5 diff
  | (name, prob) <- [ stiffish ]
  ]
  

eventTests opts solver = testGroup "Events"
  [ testCase "Exponential" $ do
      Right (eventInfo -> events) <- runKatipT ?log_env $ solver opts exponential
      length events @?= 1
      checkDiscrepancy 1e-4 (abs (eventTime (events!!0) - log 1.1))
      rootDirection (events!!0) @?= Upwards
      eventIndex (events!!0) @?= 0
  , testCase "Robertson" $ do
      let upd _ _ = vector [1.0, 0.0, 0.0]
      Right (eventInfo -> events) <- runKatipT ?log_env $ solver opts
        (snd robertson)
          { odeEvents = 
            [ EventSpec { eventCondition = \_t y -> y ! 0 - 0.0001
                        , eventUpdate = upd
                        , eventDirection = AnyDirection
                        , eventStopSolver = False
                        }
            , EventSpec { eventCondition = \_t y -> y ! 2 - 0.01
                        , eventUpdate = upd
                        , eventDirection = AnyDirection
                        , eventStopSolver = False
                        }
            ]
          , odeMaxEvents = 100
          , odeSolTimes = [0,100]
          }
      length events @?= 100
  , testCase "Bounded sine" $ do
      Right (eventInfo -> events) <- runKatipT ?log_env $ solver opts boundedSine
      length events @?= 3
      map rootDirection events @?= [Upwards, Downwards, Upwards]
      map eventIndex events @?= [0, 1, 0]
      forM_ (zip (map eventTime events) [1.119766,3.359295,5.598820]) $ \(et_got, et_exp) ->
        checkDiscrepancy 1e-4 (abs (et_exp - et_got))
  ]

{-
solveARK
  :: Katip m
  => ODEOpts CV.ODEMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveARK opts OdeProblem{..} =
  ARK.odeSolveWithEvents opts odeEvents odeMaxEvents odeRhs odeJacobian odeInitCond odeSolTimes
-}

----------------------------------------------------------------------
--                           ODE problems
----------------------------------------------------------------------

brusselator :: (String, OdeProblem)
brusselator = (,) "brusselator" $ OdeProblem
  { odeRhs = OdeRhsHaskell $ \_t x ->
      let
        u = x V.! 0
        v = x V.! 1
        w = x V.! 2
      in
      [ a - (w + 1) * u + v * u * u
      , w * u - v * u * u
      , (b - w) / eps - w * u
      ]
  , odeJacobian = Just $ \(_t :: Double) x ->
      let
        u = x V.! 0
        v = x V.! 1
        w = x V.! 2
      in (3><3)
      [ (-(w + 1.0)) + 2.0 * u * v, w - 2.0 * u * v, (-w)
      , u * u                     , (-(u * u))     , 0.0
      , (-u)                      , u              , (-1.0) / eps - u
      ]
  , odeEvents = mempty
  , odeMaxEvents = 0
  , odeInitCond = [1.2, 3.1, 3.0]
  , odeSolTimes = [0.0, 0.1 .. 10.0]
  , odeTolerances = defaultTolerances
  }
  where
    a = 1.0
    b = 3.5
    eps :: Fractional a => a
    eps = 5.0e-6

exponential = OdeProblem
  { odeRhs = OdeRhsHaskell $ \_ y -> [y V.! 0]
  , odeJacobian = Nothing
  , odeInitCond = vector [1]
  , odeEvents = events
  , odeMaxEvents = 100
  , odeSolTimes = vector [ fromIntegral k / 100 | k <- [0..(22::Int)]]
  , odeTolerances = defaultTolerances
  }
  where
    events =
      [ EventSpec { eventCondition = \_ y -> y ! 0 - 1.1
                  , eventUpdate = \_ _ -> vector [ 2 ]
                  , eventDirection = Upwards
                  , eventStopSolver = False
                  }
      ]

robertson = (,) "Robertson" $ OdeProblem
  { odeRhs = OdeRhsHaskell $ \_ (V.toList -> [y1,y2,y3]) ->
      [ -0.04 * y1 + 1.0e4 * y2 * y3
      , 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * (y2)^(2 :: Int)
      , 3.0e7 * (y2)^(2 :: Int)
      ]
  , odeJacobian = Just $ \_t (V.toList -> [_, y2, y3]) -> (3 >< 3)
      [ -0.04, 1.0e4 * y3, 1.0e4 * y2
      , 0.04, -1.0e4*y3 - 3.0e7*2*y2, -1.0e4*y2
      , 0, 3.0e7*2*y2, 0
      ]
  , odeInitCond = [1.0, 0.0, 0.0]
  , odeEvents = []
  , odeMaxEvents = 0
  , odeSolTimes = [0,5]
  , odeTolerances = defaultTolerances -- FIXME how to make this integrate indefinitely, as in the sundials example?
  }

empty = (,) "Empty system" $ OdeProblem
  { odeRhs = OdeRhsHaskell $ \_ _ -> []
  , odeJacobian = Nothing
  , odeInitCond = []
  , odeEvents = []
  , odeMaxEvents = 0
  , odeSolTimes = [0,1]
  , odeTolerances = defaultTolerances
  }

stiffish = (,) "Stiffish" $ OdeProblem
  { odeRhs = OdeRhsHaskell $ \t ((V.! 0) -> u) -> [ lamda * u + 1.0 / (1.0 + t * t) - lamda * atan t ]
  , odeJacobian = Nothing
  , odeInitCond = [0.0]
  , odeEvents = []
  , odeMaxEvents = 0
  , odeSolTimes = [0.0, 0.1 .. 10.0]
  , odeTolerances = defaultTolerances
  }
  where
    lamda = -100.0

-- A sine wave that only changes direction once it reaches ±0.9.
-- Illustrates event-specific reset function
boundedSine = OdeProblem
  { odeRhs = OdeRhsHaskell $ \_t y -> [y V.! 1, - y V.! 0]
  , odeJacobian = Nothing
  , odeInitCond = [0,1]
  , odeEvents = events
  , odeMaxEvents = 100
  , odeSolTimes = V.fromList [ 2 * pi * k / 360 | k <- [0..360]]
  , odeTolerances = defaultTolerances
  }
  where
    events =
      [ EventSpec { eventCondition = \_t y -> y ! 0 - 0.9
                     , eventUpdate = \_ y -> vector [ y ! 0, - abs (y ! 1) ]
                     , eventDirection = Upwards
                     , eventStopSolver = False
                     }
      , EventSpec { eventCondition = \_t y -> y ! 0 + 0.9
                     , eventUpdate = \_ y -> vector [ y ! 0, abs (y ! 1) ]
                     , eventDirection = Downwards
                     , eventStopSolver = False
                     }
      ]

largeTs :: V.Vector Double
largeTs = V.fromList $ 0.0 : take 12 (iterate (*10) 0.04)
