module Numeric.Sundials where

import Numeric.Sundials.Types
import Numeric.Sundials.Common
import qualified Numeric.Sundials.CVode.ODE as CV
import Control.Monad.IO.Class
import Katip

solveCV
  :: Katip m
  => ODEOpts CV.ODEMethod
  -> OdeProblem
  -> m (Either ErrorDiagnostics SundialsSolution)
solveCV = solveCommon CV.cvOdeC
