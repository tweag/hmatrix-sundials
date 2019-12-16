-- | Common infrastructure for CVode/ARKode
{-# LANGUAGE RecordWildCards #-}
module Numeric.Sundials.Common where

import Foreign.C.Types
import Numeric.Sundials.Types
import qualified Data.Vector as V
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
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
  c_local_error <- V.unsafeFreeze c_local_error
  c_var_weight <- V.unsafeFreeze c_var_weight
  c_local_error_set <- V.unsafeFreeze c_local_error_set
  return CVars {..}
