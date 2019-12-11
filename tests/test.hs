{-# LANGUAGE RecordWildCards, ScopedTypeVariables, OverloadedStrings,
             ViewPatterns, ImplicitParams #-}
import Test.Tasty
import Test.Tasty.HUnit

import qualified Numeric.Sundials.ARKode.ODE as ARK
import qualified Numeric.Sundials.CVode.ODE  as CV
import Numeric.Sundials.Types
import Numeric.LinearAlgebra as L
import qualified Data.Vector.Storable as V
import Katip
import Foreign.C.Types
import Data.Coerce
import System.IO

data OdeProblem = OdeProblem
  { odeEvents :: [EventSpec]
  , odeMaxEvents :: !Int
  , odeRhs :: OdeRhs
  , odeJacobian :: Maybe (Double -> Vector Double -> Matrix Double)
  , odeInitCond :: V.Vector Double
  , odeSolTimes :: V.Vector Double
  }

solveCV
  :: Katip m
  => ODEOpts CV.ODEMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCV opts OdeProblem{..} =
  CV.odeSolveWithEvents opts odeEvents odeMaxEvents odeRhs odeJacobian odeInitCond odeSolTimes

main = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  log_env <- registerScribe "stderr" handleScribe defaultScribeSettings =<<
    initLogEnv "test" "devel"
  let ?log_env = log_env

  defaultMain $ testGroup "Tests"
    [
      testGroup solver_name
      [ testGroup (show method) $
          let opts = defaultOpts method in
          [ withVsWithoutJacobian opts solver
          , eventTests opts solver
          ]
      | method <- [CV.BDF, CV.ADAMS]
      ]
    | (solver, solver_name) <-
      [ (solveCV, "CVode") ]
    ]

withVsWithoutJacobian opts solver = testGroup "With vs without jacobian"
  [ testCase name $ do
      Right (solutionMatrix -> solJac)   <- runKatipT ?log_env $ solver opts prob
      Right (solutionMatrix -> solNoJac) <- runKatipT ?log_env $ solver opts prob { odeJacobian = Nothing }
      assertBool "Difference too large" $ norm_2 (solJac - solNoJac) < 1e-3
  | (name, prob) <- [ brusselator ]
  ]

eventTests opts solver = testGroup "Events"
  [ testCase "Exponential" $ do
      Right (eventInfo -> events) <- runKatipT ?log_env $ solver opts exponential
      length events @?= 1
      assertBool "Difference too large" (abs (eventTime (events!!0) - log 1.1) < 1e-4)
      rootDirection (events!!0) @?= Upwards
      eventIndex (events!!0) @?= 0
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

defaultOpts :: method -> ODEOpts method
defaultOpts method = ODEOpts
  { maxNumSteps = 10000
  , minStep     = 1.0e-12
  , maxFail     = 10
  , odeMethod   = method
  , stepControl = CV.XX' 1.0e-6 1.0e-10 1 1
  , initStep    = Nothing
  }

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
      in V.fromList
      [ a - (w + 1) * u + v * u * u
      , w * u - v * u * u
      , (b - w) / eps - w * u
      ]
  , odeJacobian = Just . coerce $ \(_t :: Double) x ->
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
  , odeInitCond = V.fromList [1.2, 3.1, 3.0]
  , odeSolTimes = V.fromList [0.0, 0.1 .. 10.0]
  }
  where
    a = 1.0
    b = 3.5
    eps = 5.0e-6

exponential = OdeProblem
  { odeRhs = OdeRhsHaskell . coerce $ \(_ :: Double) y -> vector [y ! 0]
  , odeJacobian = Nothing
  , odeInitCond = vector [1]
  , odeEvents = events
  , odeMaxEvents = 100
  , odeSolTimes = vector [ fromIntegral k / 100 | k <- [0..(22::Int)]]
  }
  where
    events =
      [ EventSpec { eventCondition = \_ y -> y ! 0 - 1.1
                     , eventUpdate = \_ _ -> vector [ 2 ]
                     , eventDirection = Upwards
                     , eventStopSolver = False
                     }
      ]
