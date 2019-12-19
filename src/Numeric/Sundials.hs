module Numeric.Sundials where

import Numeric.Sundials.Types
import Numeric.Sundials.Common
import qualified Numeric.Sundials.CVode.ODE as CV
import qualified Numeric.Sundials.ARKode.ODE as ARK
import Control.Monad.IO.Class
import Katip

solveCV
  :: Katip m
  => ODEOpts CV.ODEMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCV = solveCommon CV.cvOdeC

solveARK
  :: Katip m
  => ODEOpts ARK.ODEMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveARK = solveCommon ARK.arkOdeC
