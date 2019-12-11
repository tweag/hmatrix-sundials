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
    [ testGroup (show method) . return $
        withVsWithoutJacobian (defaultOpts method) solveCV
    | method <- [CV.BDF, CV.ADAMS]
    ]

withVsWithoutJacobian opts solver = testGroup "With vs without jacobian"
  [ testCase name $ do
      Right (solutionMatrix -> solJac)   <- runKatipT ?log_env $ solver opts prob
      Right (solutionMatrix -> solNoJac) <- runKatipT ?log_env $ solver opts prob { odeJacobian = Nothing }
      assertBool "Difference too large" $ norm_2 (solJac - solNoJac) < 1e-3
  | (name, prob) <- [ brusselator ]
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

{-
brusselator :: Double -> [Double] -> [Double]
brusselator _t x = [ a - (w + 1) * u + v * u * u
                   , w * u - v * u * u
                   , (b - w) / eps - w * u
                   ]
  where
    a = 1.0
    b = 3.5
    eps = 5.0e-6
    u = x !! 0
    v = x !! 1
    w = x !! 2

brussJac :: Double -> Vector Double -> Matrix Double
brussJac _t x = tr $
  (3><3) [ (-(w + 1.0)) + 2.0 * u * v, w - 2.0 * u * v, (-w)
         , u * u                     , (-(u * u))     , 0.0
         , (-u)                      , u              , (-1.0) / eps - u
         ]
  where
    y = toList x
    u = y !! 0
    v = y !! 1
    w = y !! 2
    eps = 5.0e-6

brusselatorWithJacobian :: (MonadIO m, Katip m) => Vector Double -> Bool -> m CV.SolverResult
brusselatorWithJacobian ts usejac = CV.odeSolveRootVWith' opts
                      (OdeRhsHaskell . coerce $ \t v -> vector $ brusselator t (toList v))
                      (if usejac then Just brussJac else Nothing)
                      (vector [1.2, 3.1, 3.0])
                      [] 0
                      ts
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , maxFail     = 10
                   , odeMethod = CV.BDF
                   , stepControl = CV.XX' 1.0e-6 1.0e-10 1 1
                   , initStep = Nothing
                   }

-}
