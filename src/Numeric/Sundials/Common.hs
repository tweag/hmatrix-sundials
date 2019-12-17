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
import Numeric.LinearAlgebra.HMatrix hiding (Vector)
import GHC.Prim
import Control.Monad.IO.Class
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
  c_diagnostics <- VSM.new 11
  c_root_info <- VSM.new $ V.length odeEvents
  c_event_index <- VSM.new odeMaxEvents
  c_event_time <- VSM.new odeMaxEvents
  c_actual_event_direction <- VSM.new odeMaxEvents
  c_n_events <- VSM.new 1
  c_n_rows <- VSM.new 1
  c_local_error <- VSM.new $ VS.length odeInitCond
  c_var_weight <- VSM.new $ VS.length odeInitCond
  c_local_error_set <- VSM.new 1
  return CVars {..}

-- NB: the mutable CVars must not be used after this
freezeCVars :: CVars (V.MVector RealWorld) -> IO (CVars V.Vector)
freezeCVars CVars{..} = do
  c_diagnostics <- V.unsafeFreeze c_diagnostics
  c_root_info <- V.unsafeFreeze c_root_info
  c_event_index <- V.unsafeFreeze c_event_index
  c_event_time <- V.unsafeFreeze c_event_time
  c_actual_event_direction <- V.unsafeFreeze c_actual_event_direction
  c_n_events <- V.unsafeFreeze c_n_events
  c_n_rows <- V.unsafeFreeze c_n_rows
  c_output_mat <- V.unsafeFreeze c_output_mat
  c_local_error <- V.unsafeFreeze c_local_error
  c_var_weight <- V.unsafeFreeze c_var_weight
  c_local_error_set <- V.unsafeFreeze c_local_error_set
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

-- | The common solving logic between ARKode and CVode
solveCommon
  :: Katip m
  => ODEOpts method
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCommon ODEOpts{..} OdeProblem{..}

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
  (rhs_funptr :: FunPtr OdeRhsCType, userdata :: Ptr UserData) <-
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
        funptr <- mkOdeRhsC funIO
        return (funptr, nullPtr)
  undefined

foreign import ccall "wrapper"
  mkOdeRhsC :: OdeRhsCType -> IO (FunPtr OdeRhsCType)
