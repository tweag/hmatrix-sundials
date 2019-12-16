-- | Common infrastructure for CVode/ARKode
{-# LANGUAGE RecordWildCards #-}
module Numeric.Sundials.Common where

import Foreign.C.Types
import Numeric.Sundials.Types
import qualified Data.Vector as VB
import qualified Data.Vector.Storable as V
import qualified Data.Vector.Storable.Mutable as VM
import GHC.Prim

-- | A collection of variables that we allocate on the Haskell side and
-- pass into the C code to be filled.
data CVars vec = CVars
  { c_diagnostics :: vec CInt
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

allocateCVars :: OdeProblem -> IO (CVars (V.MVector RealWorld))
allocateCVars OdeProblem{..} = do 
  c_diagnostics <- VM.new 11
  c_root_info <- VM.new $ VB.length odeEvents
  c_event_index <- VM.new odeMaxEvents
  c_event_time <- VM.new odeMaxEvents
  c_actual_event_direction <- VM.new odeMaxEvents
  c_n_events <- VM.new 1
  c_n_rows <- VM.new 1
  c_local_error <- VM.new $ V.length odeInitCond
  c_var_weight <- VM.new $ V.length odeInitCond
  c_local_error_set <- VM.new 1
  return CVars {..}
