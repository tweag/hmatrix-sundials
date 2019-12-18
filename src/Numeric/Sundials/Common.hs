-- | Common infrastructure for CVode/ARKode
{-# LANGUAGE TemplateHaskell #-}
module Numeric.Sundials.Common where

import Foreign.C.Types
import Foreign.Ptr
import Foreign.Storable (peek, poke)
import Numeric.Sundials.Types
import qualified Numeric.Sundials.Foreign as T
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Data.Maybe
import Numeric.LinearAlgebra.HMatrix hiding (Vector)
import GHC.Prim
import Control.Monad.IO.Class
import Control.Monad.Cont
import Control.Exception
import Katip
import Language.Haskell.TH

-- | A collection of variables that we allocate on the Haskell side and
-- pass into the C code to be filled.
data CVars vec = CVars
  { c_diagnostics :: vec SunIndexType
    -- ^ Mutable vector to which we write diagnostic data while
    -- solving. Its size corresponds to the number of fields in
    -- 'SundialsDiagnostics'.
  , c_root_info :: vec CInt
    -- ^ Just a temporary vector (of the size equal to the number of event
    -- specs) that we use to get root info. Isn't used for output.
  , c_event_index :: vec CInt
    -- ^ For each event occurrence, this indicates which of the events
    -- occurred. Size: max_num_events.
  , c_event_time :: vec CDouble
    -- ^ For each event occurrence, this indicates the time of the
    -- occurrence. Size: max_num_events.
  , c_n_events :: vec CInt
    -- ^ Vector of size 1 that gives the total number of events occurred.
  , c_n_rows :: vec CInt
    -- ^ The total number of rows in the output matrix.
  , c_output_mat :: vec CDouble
    -- ^ The output matrix stored in the row-major order.
    -- Dimensions: (1 + dim) * (2 * max_events + nTs).
  , c_actual_event_direction :: vec CInt
    -- ^ Vector of size max_num_events that gives the direction of the
    -- occurred event.
  , c_local_error :: vec CDouble
    -- ^ Vector containing local error estimates. Size: the dimensionality
    -- of the system.
  , c_var_weight :: vec CDouble
    -- ^ Vector containing variable weights (derived from the tolerances).
    -- Size: the dimensionality of the system.
  , c_local_error_set :: vec CInt
    -- The flag (size 1) indicating whether c_local_error is filled with meaningful
    -- values. *Should be initialized with 0.*
  }

allocateCVars :: OdeProblem -> IO (CVars (VS.MVector RealWorld))
allocateCVars OdeProblem{..} = do 
  let dim = VS.length odeInitCond
  c_diagnostics <- VSM.new 11
  c_root_info <- VSM.new $ V.length odeEvents
  c_event_index <- VSM.new odeMaxEvents
  c_event_time <- VSM.new odeMaxEvents
  c_actual_event_direction <- VSM.new odeMaxEvents
  c_n_events <- VSM.new 1
  c_n_rows <- VSM.new 1
  c_local_error <- VSM.new dim
  c_var_weight <- VSM.new dim
  c_local_error_set <- VSM.new 1
  c_output_mat <- VSM.new $
    (1 + dim) * (2 * odeMaxEvents + VS.length odeSolTimes)
  return CVars {..}

-- NB: the mutable CVars must not be used after this
freezeCVars :: CVars (VS.MVector RealWorld) -> IO (CVars VS.Vector)
freezeCVars CVars{..} = do
  c_diagnostics <- VS.unsafeFreeze c_diagnostics
  c_root_info <- VS.unsafeFreeze c_root_info
  c_event_index <- VS.unsafeFreeze c_event_index
  c_event_time <- VS.unsafeFreeze c_event_time
  c_actual_event_direction <- VS.unsafeFreeze c_actual_event_direction
  c_n_events <- VS.unsafeFreeze c_n_events
  c_n_rows <- VS.unsafeFreeze c_n_rows
  c_output_mat <- VS.unsafeFreeze c_output_mat
  c_local_error <- VS.unsafeFreeze c_local_error
  c_var_weight <- VS.unsafeFreeze c_var_weight
  c_local_error_set <- VS.unsafeFreeze c_local_error_set
  return CVars {..}

-- | Similar to 'CVars', except these are immutable values that are
-- accessed (read-only) by the C code and specify the system to be solved.
data CConsts = CConsts
  { c_dim :: SunIndexType -- ^ the dimensionality (number of variables/equations)
  , c_method :: CInt -- ^ the ODE method (specific to the solver)
  , c_n_sol_times :: CInt
  , c_sol_time :: VS.Vector CDouble
  , c_init_cond :: VS.Vector CDouble
  , c_rhs :: FunPtr OdeRhsCType
  , c_rhs_userdata :: Ptr UserData
  , c_rtol :: CDouble
  , c_atol :: VS.Vector CDouble
  , c_n_event_specs :: CInt
  , c_event_fn :: CDouble -> Ptr T.SunVector -> Ptr CDouble -> Ptr () -> IO CInt
  , c_apply_event :: CInt -> CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> IO CInt
  , c_jac_set :: CInt
  , c_jac :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunMatrix
          -> Ptr () -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunVector
          -> IO CInt
  , c_requested_event_direction :: VS.Vector CInt
  , c_event_stops_solver :: VS.Vector CInt
  , c_max_events :: CInt
  , c_minstep :: CDouble
  , c_max_n_steps :: SunIndexType
  , c_max_err_test_fails :: CInt
  , c_init_step_size_set :: CInt
  , c_init_step_size :: CDouble
  }

class Method method where
  methodToInt :: method -> CInt

withCConsts
  :: Method method
  => ODEOpts method
  -> OdeProblem
  -> (CConsts -> IO r)
  -> IO r
withCConsts ODEOpts{..} OdeProblem{..} = runContT $ do
  let
    dim = VS.length c_init_cond
    c_init_cond = coerce odeInitCond
    c_dim = fromIntegral c_dim
    c_n_sol_times = fromIntegral . VS.length $ odeSolTimes
    c_sol_time = coerce odeSolTimes
    c_rtol = relTolerance odeTolerances
    c_atol = either (VS.replicate dim) id $ absTolerances odeTolerances
    c_minstep = coerce minStep
    c_max_n_steps = fromIntegral maxNumSteps
    c_max_err_test_fails = fromIntegral maxFail
    c_init_step_size_set = fromIntegral . fromEnum $ isJust initStep
    c_init_step_size = coerce . fromMaybe undefined $ initStep
    c_n_event_specs = fromIntegral $ V.length odeEvents
    c_requested_event_direction = V.convert $ V.map (directionToInt . eventDirection) odeEvents
    c_event_fn t y_ptr out_ptr _ptr = do
      y <- sunVecVals <$> peek y_ptr
      let vals = V.convert $ V.map (\ev -> coerce (eventCondition ev) t y) odeEvents
      -- FIXME: We should be able to use poke somehow
      T.vectorToC vals (fromIntegral c_n_event_specs) out_ptr
      return 0
    c_apply_event event_index t y_ptr y'_ptr = do
      y_vec <- peek y_ptr
      let
        ev = odeEvents V.! (fromIntegral event_index)
        y' = coerce (eventUpdate ev) t (sunVecVals y_vec)
      poke y'_ptr $ SunVector
        { sunVecN = sunVecN y_vec
        , sunVecVals = y'
        }
      return 0
    c_event_stops_solver = 
      V.convert
      . V.map (fromIntegral . fromEnum . eventStopSolver)
      $ odeEvents
    c_max_events = fromIntegral odeMaxEvents
    c_jac_set = fromIntegral . fromEnum $ isJust odeJacobian
    c_jac t y _fy jacS _ptr _tmp1 _tmp2 _tmp3 = do
      case odeJacobian of
        Nothing   -> undefined
        Just jacI -> do j <- matrixToSunMatrix . jacI (coerce t) <$> (coerce $ sunVecVals <$> peek y)
                        poke jacS j
                        return 0
    c_method = methodToInt odeMethod

  (c_rhs, c_rhs_userdata) <-
    case odeRhs of
      OdeRhsC ptr u -> return (ptr, u)
      OdeRhsHaskell fun -> do
        let
          funIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr UserData -> IO CInt
          funIO t y f _ptr = do
            sv <- peek y
            poke f $ SunVector { sunVecN = sunVecN sv
                               , sunVecVals = fun t (sunVecVals sv)
                               }
            return 0
        funptr <- ContT $ bracket (mkOdeRhsC funIO) freeHaskellFunPtr
        return (funptr, nullPtr)
  return CConsts{..}

matrixToSunMatrix :: Matrix Double -> T.SunMatrix
matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
  where
    nr = fromIntegral $ rows m
    nc = fromIntegral $ cols m
    -- FIXME: efficiency
    vs = VS.fromList $ map coerce $ concat $ toLists m

-- Contrary to the documentation, it appears that CVodeGetRootInfo
-- may use both 1 and -1 to indicate a root, depending on the
-- direction of the sign change. See near the end of cvRootfind.
intToDirection :: Integral d => d -> Maybe CrossingDirection
intToDirection d =
  case d of
    1  -> Just Upwards
    -1 -> Just Downwards
    _  -> Nothing

-- | Almost inverse of 'intToDirection'. Map 'Upwards' to 1, 'Downwards' to
-- -1, and 'AnyDirection' to 0.
directionToInt :: Integral d => CrossingDirection -> d
directionToInt d =
  case d of
    Upwards -> 1
    Downwards -> -1
    AnyDirection -> 0

foreign import ccall "wrapper"
  mkOdeRhsC :: OdeRhsCType -> IO (FunPtr OdeRhsCType)

assembleSolverResult
  :: CInt
  -> CVars VS.Vector
  -> IO (Either ErrorDiagnostics SundialsSolution)
assembleSolverResult = undefined

-- | The common solving logic between ARKode and CVode
solveCommon
  :: (Method method, Katip m)
  => (CConsts -> CVars (VS.MVector RealWorld) -> ReportErrorFn -> IO CInt)
      -- ^ the CVode/ARKode solving function; mostly inline-C code
  -> ODEOpts method
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCommon solve_c opts problem@(OdeProblem{..})

  | VS.null odeInitCond = -- 0-dimensional (empty) system

    return . Right $ SundialsSolution
      { actualTimeGrid = odeSolTimes
      , solutionMatrix = (VS.length odeSolTimes >< 0) []
      , eventInfo = []
      , diagnostics = emptyDiagnostics
      }

  | otherwise = do

    report_error <- logWithKatip
    liftIO $ do -- the rest is in the IO monad
    vars <- allocateCVars problem
    ret <- withCConsts opts problem $ \consts ->
      solve_c consts vars report_error
    frozenVars <- freezeCVars vars
    assembleSolverResult ret frozenVars
