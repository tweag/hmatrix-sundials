{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FlexibleContexts #-}

-- |
-- Solution of ordinary differential equation (ODE) initial value problems.
-- See <https://computation.llnl.gov/projects/sundials/sundials-software> for more detail.
module Numeric.Sundials.ARKode.ODE
  ( odeSolveWithEvents
  , ODEMethod(..)
  , StepControl(..)
  ) where

import qualified Language.C.Inline as C
import qualified Language.C.Inline.Unsafe as CU

import           Data.Monoid ((<>))
import           Data.Maybe (isJust)

import           Foreign.C.Types (CDouble, CInt)
import           Foreign.Ptr (Ptr)
import           Foreign.Storable (poke, peek)

import qualified Data.Vector.Storable as V

import           Data.Coerce (coerce)
import           GHC.Generics (C1, Constructor, (:+:)(..), D1, Rep, Generic, M1(..),
                               from, conName)

import           Numeric.LinearAlgebra.Devel (createVector)

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix, rows,
                                                cols, toLists, size, reshape,
                                                (><))

import           Numeric.Sundials.Types
import qualified Numeric.Sundials.Arkode as T
import           Numeric.Sundials.Arkode (sDIRK_2_1_2,
                                          bILLINGTON_3_3_2,
                                          tRBDF2_3_3_2,
                                          kVAERNO_4_2_3,
                                          aRK324L2SA_DIRK_4_2_3,
                                          cASH_5_2_4,
                                          cASH_5_3_4,
                                          sDIRK_5_3_4,
                                          kVAERNO_5_3_4,
                                          aRK436L2SA_DIRK_6_3_4,
                                          kVAERNO_7_4_5,
                                          aRK548L2SA_DIRK_8_4_5,
                                          hEUN_EULER_2_1_2,
                                          bOGACKI_SHAMPINE_4_2_3,
                                          aRK324L2SA_ERK_4_2_3,
                                          zONNEVELD_5_3_4,
                                          aRK436L2SA_ERK_6_3_4,
                                          sAYFY_ABURUB_6_3_4,
                                          cASH_KARP_6_4_5,
                                          fEHLBERG_6_4_5,
                                          dORMAND_PRINCE_7_4_5,
                                          aRK548L2SA_ERK_8_4_5,
                                          vERNER_8_5_6,
                                          fEHLBERG_13_7_8)

import Katip
import Control.Monad.IO.Class


C.context (C.baseCtx <> C.vecCtx <> C.funCtx <> sunCtx)

C.include "<stdlib.h>"
C.include "<stdio.h>"
C.include "<math.h>"
C.include "<arkode/arkode.h>"
C.include "<arkode/arkode_ls.h>"
C.include "<arkode/arkode_arkstep.h>"
C.include "<nvector/nvector_serial.h>"
C.include "<sunmatrix/sunmatrix_dense.h>"
C.include "<sunlinsol/sunlinsol_dense.h>"
C.include "<sundials/sundials_types.h>"
C.include "<sundials/sundials_math.h>"
C.include "Numeric/Sundials/Arkode_hsc.h"
C.include "../../../helpers.h"


-- | Stepping functions
data ODEMethod = SDIRK_2_1_2            Jacobian
               | SDIRK_2_1_2'
               | BILLINGTON_3_3_2       Jacobian
               | BILLINGTON_3_3_2'
               | TRBDF2_3_3_2           Jacobian
               | TRBDF2_3_3_2'
               | KVAERNO_4_2_3          Jacobian
               | KVAERNO_4_2_3'
               | ARK324L2SA_DIRK_4_2_3  Jacobian
               | ARK324L2SA_DIRK_4_2_3'
               | CASH_5_2_4             Jacobian
               | CASH_5_2_4'
               | CASH_5_3_4             Jacobian
               | CASH_5_3_4'
               | SDIRK_5_3_4            Jacobian
               | SDIRK_5_3_4'
               | KVAERNO_5_3_4          Jacobian
               | KVAERNO_5_3_4'
               | ARK436L2SA_DIRK_6_3_4  Jacobian
               | ARK436L2SA_DIRK_6_3_4'
               | KVAERNO_7_4_5          Jacobian
               | KVAERNO_7_4_5'
               | ARK548L2SA_DIRK_8_4_5  Jacobian
               | ARK548L2SA_DIRK_8_4_5'
               | HEUN_EULER_2_1_2         Jacobian
               | HEUN_EULER_2_1_2'
               | BOGACKI_SHAMPINE_4_2_3   Jacobian
               | BOGACKI_SHAMPINE_4_2_3'
               | ARK324L2SA_ERK_4_2_3     Jacobian
               | ARK324L2SA_ERK_4_2_3'
               | ZONNEVELD_5_3_4          Jacobian
               | ZONNEVELD_5_3_4'
               | ARK436L2SA_ERK_6_3_4     Jacobian
               | ARK436L2SA_ERK_6_3_4'
               | SAYFY_ABURUB_6_3_4       Jacobian
               | SAYFY_ABURUB_6_3_4'
               | CASH_KARP_6_4_5          Jacobian
               | CASH_KARP_6_4_5'
               | FEHLBERG_6_4_5         Jacobian
               | FEHLBERG_6_4_5'
               | DORMAND_PRINCE_7_4_5     Jacobian
               | DORMAND_PRINCE_7_4_5'
               | ARK548L2SA_ERK_8_4_5     Jacobian
               | ARK548L2SA_ERK_8_4_5'
               | VERNER_8_5_6            Jacobian
               | VERNER_8_5_6'
               | FEHLBERG_13_7_8         Jacobian
               | FEHLBERG_13_7_8'
  deriving Generic

constrName :: (HasConstructor (Rep a), Generic a)=> a -> String
constrName = genericConstrName . from

class HasConstructor (f :: * -> *) where
  genericConstrName :: f x -> String

instance HasConstructor f => HasConstructor (D1 c f) where
  genericConstrName (M1 x) = genericConstrName x

instance (HasConstructor x, HasConstructor y) => HasConstructor (x :+: y) where
  genericConstrName (L1 l) = genericConstrName l
  genericConstrName (R1 r) = genericConstrName r

instance Constructor c => HasConstructor (C1 c f) where
  genericConstrName x = conName x

instance Show ODEMethod where
  show x = constrName x

-- FIXME: We can probably do better here with generics
getMethod :: ODEMethod -> Int
getMethod (SDIRK_2_1_2 _)            = sDIRK_2_1_2
getMethod (SDIRK_2_1_2')             = sDIRK_2_1_2
getMethod (BILLINGTON_3_3_2 _)       = bILLINGTON_3_3_2
getMethod (BILLINGTON_3_3_2')        = bILLINGTON_3_3_2
getMethod (TRBDF2_3_3_2 _)           = tRBDF2_3_3_2
getMethod (TRBDF2_3_3_2')            = tRBDF2_3_3_2
getMethod (KVAERNO_4_2_3  _)         = kVAERNO_4_2_3
getMethod (KVAERNO_4_2_3')           = kVAERNO_4_2_3
getMethod (ARK324L2SA_DIRK_4_2_3 _)  = aRK324L2SA_DIRK_4_2_3
getMethod (ARK324L2SA_DIRK_4_2_3')   = aRK324L2SA_DIRK_4_2_3
getMethod (CASH_5_2_4 _)             = cASH_5_2_4
getMethod (CASH_5_2_4')              = cASH_5_2_4
getMethod (CASH_5_3_4 _)             = cASH_5_3_4
getMethod (CASH_5_3_4')              = cASH_5_3_4
getMethod (SDIRK_5_3_4 _)            = sDIRK_5_3_4
getMethod (SDIRK_5_3_4')             = sDIRK_5_3_4
getMethod (KVAERNO_5_3_4 _)          = kVAERNO_5_3_4
getMethod (KVAERNO_5_3_4')           = kVAERNO_5_3_4
getMethod (ARK436L2SA_DIRK_6_3_4 _)  = aRK436L2SA_DIRK_6_3_4
getMethod (ARK436L2SA_DIRK_6_3_4')   = aRK436L2SA_DIRK_6_3_4
getMethod (KVAERNO_7_4_5 _)          = kVAERNO_7_4_5
getMethod (KVAERNO_7_4_5')           = kVAERNO_7_4_5
getMethod (ARK548L2SA_DIRK_8_4_5 _)  = aRK548L2SA_DIRK_8_4_5
getMethod (ARK548L2SA_DIRK_8_4_5')   = aRK548L2SA_DIRK_8_4_5
getMethod (HEUN_EULER_2_1_2 _)       = hEUN_EULER_2_1_2
getMethod (HEUN_EULER_2_1_2')        = hEUN_EULER_2_1_2
getMethod (BOGACKI_SHAMPINE_4_2_3 _) = bOGACKI_SHAMPINE_4_2_3
getMethod (BOGACKI_SHAMPINE_4_2_3')  = bOGACKI_SHAMPINE_4_2_3
getMethod (ARK324L2SA_ERK_4_2_3 _)   = aRK324L2SA_ERK_4_2_3
getMethod (ARK324L2SA_ERK_4_2_3')    = aRK324L2SA_ERK_4_2_3
getMethod (ZONNEVELD_5_3_4 _)        = zONNEVELD_5_3_4
getMethod (ZONNEVELD_5_3_4')         = zONNEVELD_5_3_4
getMethod (ARK436L2SA_ERK_6_3_4 _)   = aRK436L2SA_ERK_6_3_4
getMethod (ARK436L2SA_ERK_6_3_4')    = aRK436L2SA_ERK_6_3_4
getMethod (SAYFY_ABURUB_6_3_4 _)     = sAYFY_ABURUB_6_3_4
getMethod (SAYFY_ABURUB_6_3_4')      = sAYFY_ABURUB_6_3_4
getMethod (CASH_KARP_6_4_5 _)        = cASH_KARP_6_4_5
getMethod (CASH_KARP_6_4_5')         = cASH_KARP_6_4_5
getMethod (FEHLBERG_6_4_5 _)         = fEHLBERG_6_4_5
getMethod (FEHLBERG_6_4_5' )         = fEHLBERG_6_4_5
getMethod (DORMAND_PRINCE_7_4_5 _)   = dORMAND_PRINCE_7_4_5
getMethod (DORMAND_PRINCE_7_4_5')    = dORMAND_PRINCE_7_4_5
getMethod (ARK548L2SA_ERK_8_4_5 _)   = aRK548L2SA_ERK_8_4_5
getMethod (ARK548L2SA_ERK_8_4_5')    = aRK548L2SA_ERK_8_4_5
getMethod (VERNER_8_5_6 _)           = vERNER_8_5_6
getMethod (VERNER_8_5_6')            = vERNER_8_5_6
getMethod (FEHLBERG_13_7_8 _)        = fEHLBERG_13_7_8
getMethod (FEHLBERG_13_7_8')         = fEHLBERG_13_7_8

getJacobian :: ODEMethod -> Maybe Jacobian
getJacobian (SDIRK_2_1_2 j)            = Just j
getJacobian (BILLINGTON_3_3_2 j)       = Just j
getJacobian (TRBDF2_3_3_2 j)           = Just j
getJacobian (KVAERNO_4_2_3  j)         = Just j
getJacobian (ARK324L2SA_DIRK_4_2_3 j)  = Just j
getJacobian (CASH_5_2_4 j)             = Just j
getJacobian (CASH_5_3_4 j)             = Just j
getJacobian (SDIRK_5_3_4 j)            = Just j
getJacobian (KVAERNO_5_3_4 j)          = Just j
getJacobian (ARK436L2SA_DIRK_6_3_4 j)  = Just j
getJacobian (KVAERNO_7_4_5 j)          = Just j
getJacobian (ARK548L2SA_DIRK_8_4_5 j)  = Just j
getJacobian (HEUN_EULER_2_1_2 j)       = Just j
getJacobian (BOGACKI_SHAMPINE_4_2_3 j) = Just j
getJacobian (ARK324L2SA_ERK_4_2_3 j)   = Just j
getJacobian (ZONNEVELD_5_3_4 j)        = Just j
getJacobian (ARK436L2SA_ERK_6_3_4 j)   = Just j
getJacobian (SAYFY_ABURUB_6_3_4 j)     = Just j
getJacobian (CASH_KARP_6_4_5 j)        = Just j
getJacobian (FEHLBERG_6_4_5 j)         = Just j
getJacobian (DORMAND_PRINCE_7_4_5 j)   = Just j
getJacobian (ARK548L2SA_ERK_8_4_5 j)   = Just j
getJacobian (VERNER_8_5_6 j)           = Just j
getJacobian (FEHLBERG_13_7_8 j)        = Just j
getJacobian _                          = Nothing

odeSolveVWith' :: (MonadIO m, Katip m)
  => ODEOpts ODEMethod
  -> ODEMethod
  -> StepControl
  -> Maybe Double -- ^ initial step size - by default, ARKode
                  -- estimates the initial step size to be the
                  -- solution \(h\) of the equation
                  -- \(\|\frac{h^2\ddot{y}}{2}\| = 1\), where
                  -- \(\ddot{y}\) is an estimated value of the second
                  -- derivative of the solution at \(t_0\)
  -> (Double -> V.Vector Double -> V.Vector Double) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector Double                     -- ^ Initial conditions
  -> V.Vector Double                     -- ^ Desired solution times
  -> m (Either (Matrix Double, Int) (Matrix Double, SundialsDiagnostics)) -- ^ Error code or solution
odeSolveVWith' opts method control initStepSize f y0 tt = do
  r <- solveOdeC (fromIntegral $ maxFail opts)
                 (fromIntegral $ maxNumSteps opts) (coerce $ minStep opts)
                 (fromIntegral $ getMethod method) (coerce initStepSize) jacH (scise control)
                 (coerce f) (coerce y0) (coerce tt)
  return $ case r of
    Left  (v, c) -> Left  (reshape l (coerce v), fromIntegral c)
    Right (v, d)
      | V.null y0 -> Right ((V.length tt >< 0) [], emptyDiagnostics)
      | otherwise -> Right (reshape l (coerce v), d)
  where
    l = size y0
    scise (X aTol rTol)                          = coerce (V.replicate l aTol, rTol)
    scise (X' aTol rTol)                         = coerce (V.replicate l aTol, rTol)
    scise (XX' aTol rTol yScale _yDotScale)      = coerce (V.replicate l aTol, yScale * rTol)
    -- FIXME; Should we check that the length of ss is correct?
    scise (ScXX' aTol rTol yScale _yDotScale ss) = coerce (V.map (* aTol) ss, yScale * rTol)
    jacH = fmap (\g t v -> matrixToSunMatrix $ g (coerce t) (coerce v)) $
           getJacobian method
    matrixToSunMatrix m = T.SunMatrix { T.rows = nr, T.cols = nc, T.vals = vs }
      where
        nr = fromIntegral $ rows m
        nc = fromIntegral $ cols m
        -- FIXME: efficiency
        vs = V.fromList $ map coerce $ concat $ toLists m

-- | This function implements the same interface as
-- 'Numeric.Sundials.CVode.ODE.odeSolveWithEvents', although it does not
-- currently support events.
odeSolveWithEvents
  :: (MonadIO m, Katip m)
  => ODEOpts ODEMethod
  -> [EventSpec]
    -- ^ Event specifications
  -> Int
    -- ^ Maximum number of events
  -> (Double -> V.Vector Double -> V.Vector Double)
    -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> Maybe (Double -> Vector Double -> Matrix Double)
    -- ^ The Jacobian (optional)
  -> V.Vector Double
    -- ^ Initial conditions
  -> V.Vector Double
    -- ^ Desired solution times
  -> m (Either Int SundialsSolution)
    -- ^ Either an error code or a solution
odeSolveWithEvents opts events _ rhs _mb_jac y0 times
  | (not . null) events =
      -- Call error rather than return a Left because this is a programming
      -- error, not just a runtime issue.
      error $ "ARKode called with a non-empty list of events (" ++ show (length events) ++
      " in total).\
      \ ARKode does not support events at this point and should not be passed any."
  | otherwise = do
      result :: Either (Matrix Double, Int)
                       (Matrix Double, SundialsDiagnostics)
                       <-
        odeSolveVWith' opts
          (odeMethod opts)
          (stepControl opts)
          (initStep opts)
          rhs y0 times
      return $
        case result of
          Left (_, code) -> Left code
          Right (mx, diagn) ->
            Right $ SundialsSolution
                { actualTimeGrid = times
                , solutionMatrix =
                    -- Note: at this time, ARKode's output matrix does not
                    -- include the time column, so we're not dropping it
                    -- here unlike in CVode. If/when we add event support
                    -- to ARKode, this is going to change.
                    mx
                , eventInfo = []
                , diagnostics = diagn
                }

solveOdeC :: (MonadIO m, Katip m) =>
  CInt ->
  T.SunIndexType ->
  CDouble ->
  CInt ->
  Maybe CDouble ->
  (Maybe (CDouble -> V.Vector CDouble -> T.SunMatrix)) ->
  (V.Vector CDouble, CDouble) ->
  (CDouble -> V.Vector CDouble -> V.Vector CDouble) -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector CDouble -- ^ Initial conditions
  -> V.Vector CDouble -- ^ Desired solution times
  -> m (Either (V.Vector CDouble, CInt) (V.Vector CDouble, SundialsDiagnostics)) -- ^ Partial solution and error code or
                                                                             -- solution and diagnostics
solveOdeC maxErrTestFails maxNumSteps_ minStep_ method initStepSize
          jacH (aTols, rTol) fun f0 ts
  | V.null f0 = -- 0-dimensional (empty) system
    return $ Right (V.empty, emptyDiagnostics)
  | otherwise = do
  report_error <- logWithKatip
  liftIO $ do
  let isInitStepSize :: CInt
      isInitStepSize = fromIntegral $ fromEnum $ isJust initStepSize
      ss :: CDouble
      ss = case initStepSize of
             -- It would be better to put an error message here but
             -- inline-c seems to evaluate this even if it is never
             -- used :(
             Nothing -> 0.0
             Just x  -> x

  let dim = V.length f0
      nEq :: T.SunIndexType
      nEq = fromIntegral dim
      nTs :: CInt
      nTs = fromIntegral $ V.length ts
  quasiMatrixRes <- createVector ((fromIntegral dim) * (fromIntegral nTs))
  qMatMut <- V.thaw quasiMatrixRes
  diagn :: V.Vector T.SunIndexType <- createVector 10 -- FIXME
  diagMut <- V.thaw diagn
  -- We need the types that sundials expects. These are tied together
  -- in 'CLangToHaskellTypes'. FIXME: The Haskell type is currently empty!
  let funIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr () -> IO CInt
      funIO t y f _ptr = do
        sv <- peek y
        poke f $ T.SunVector { T.sunVecN = T.sunVecN sv
                             , T.sunVecVals = fun t (T.sunVecVals sv)
                             }
        [CU.exp| int{ 0 } |]
  let isJac :: CInt
      isJac = fromIntegral $ fromEnum $ isJust jacH
      jacIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunMatrix ->
               Ptr () -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunVector ->
               IO CInt
      jacIO t y _fy jacS _ptr _tmp1 _tmp2 _tmp3 = do
        case jacH of
          Nothing   -> error "Numeric.Sundials.ARKode.ODE: Jacobian not defined"
          Just jacI -> do j <- jacI t <$> (T.sunVecVals <$> peek y)
                          poke jacS j
                          -- FIXME: I don't understand what this comment means
                          -- Unsafe since the function will be called many times.
                          [CU.exp| int{ 0 } |]

  res <- [C.block| int {
                         /* general problem variables */

                         int flag;                  /* reusable error-checking flag                 */
                         int i, j;                  /* reusable loop indices                        */
                         N_Vector y = NULL;         /* empty vector for storing solution            */
                         N_Vector tv = NULL;        /* empty vector for storing absolute tolerances */
                         SUNMatrix A = NULL;        /* empty matrix for linear solver               */
                         SUNLinearSolver LS = NULL; /* empty linear solver object                   */
                         void *arkode_mem = NULL;   /* empty ARKode memory structure                */
                         realtype t;
                         long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;

                         /* general problem parameters */

                         realtype T0 = RCONST(($vec-ptr:(double *ts))[0]); /* initial time              */
                         sunindextype NEQ = $(sunindextype nEq);             /* number of dependent vars. */

                         /* Initialize data structures */

                         ARKErrHandlerFn report_error = $fun:(void (*report_error)(int,const char*, const char*, char*, void*));

                         y = N_VNew_Serial(NEQ); /* Create serial vector for solution */
                         if (check_flag((void *)y, "N_VNew_Serial", 0, report_error)) return 1;
                         /* Specify initial condition */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(y,i) = ($vec-ptr:(double *f0))[i];
                         };

                         tv = N_VNew_Serial(NEQ); /* Create serial vector for absolute tolerances */
                         if (check_flag((void *)tv, "N_VNew_Serial", 0, report_error)) return 1;
                         /* Specify tolerances */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(tv,i) = ($vec-ptr:(double *aTols))[i];
                         };

                         /* Call ARKStepCreate to initialize the integrator memory and specify the */
                         /* right-hand side function in y'=f(t,y), the inital time T0, and      */
                         /* the initial dependent variable vector y. */

                         /* Here we use the C types defined in helpers.h which tie up with */
                         /* the Haskell types defined in CLangToHaskellTypes                             */
                         if ($(int method) < MIN_DIRK_NUM) {
                           arkode_mem = ARKStepCreate($fun:(int (* funIO) (double t, SunVector y[], SunVector dydt[], void * params)), NULL, T0, y);
                           if (check_flag(arkode_mem, "ARKStepCreate", 0, report_error)) return 1;
                         } else {
                           arkode_mem = ARKStepCreate(NULL, $fun:(int (* funIO) (double t, SunVector y[], SunVector dydt[], void * params)), T0, y);
                           if (check_flag(arkode_mem, "ARKStepCreate", 0, report_error)) return 1;
                         }

                         flag = ARKStepSetMinStep(arkode_mem, $(double minStep_));
                         if (check_flag(&flag, "ARKStepSetMinStep", 1, report_error)) return 1;
                         flag = ARKStepSetMaxNumSteps(arkode_mem, $(sunindextype maxNumSteps_));
                         if (check_flag(&flag, "ARKStepSetMaxNumSteps", 1, report_error)) return 1;
                         flag = ARKStepSetMaxErrTestFails(arkode_mem, $(int maxErrTestFails));
                         if (check_flag(&flag, "ARKStepSetMaxErrTestFails", 1, report_error)) return 1;

                         /* Set routines */
                         flag = ARKStepSVtolerances(arkode_mem, $(double rTol), tv);
                         if (check_flag(&flag, "ARKStepSVtolerances", 1, report_error)) return 1;

                         /* Initialize dense matrix data structure and solver */
                         A = SUNDenseMatrix(NEQ, NEQ);
                         if (check_flag((void *)A, "SUNDenseMatrix", 0, report_error)) return 1;
                         LS = SUNDenseLinearSolver(y, A);
                         if (check_flag((void *)LS, "SUNDenseLinearSolver", 0, report_error)) return 1;

                         /* Attach matrix and linear solver */
                         flag = ARKStepSetLinearSolver(arkode_mem, LS, A);
                         if (check_flag(&flag, "ARKStepSetLinearSolver", 1, report_error)) return 1;

                         /* Set the initial step size if there is one */
                         if ($(int isInitStepSize)) {
                           /* FIXME: We could check if the initial step size is 0 */
                           /* or even NaN and then throw an error                 */
                           flag = ARKStepSetInitStep(arkode_mem, $(double ss));
                           if (check_flag(&flag, "ARKStepSetInitStep", 1, report_error)) return 1;
                         }

                         /* Set the Jacobian if there is one */
                         if ($(int isJac)) {
                           flag = ARKStepSetJacFn(arkode_mem, $fun:(int (* jacIO) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
                           if (check_flag(&flag, "ARKStepSetJacFn", 1, report_error)) return 1;
                         }

                         /* Store initial conditions */
                         for (j = 0; j < NEQ; j++) {
                           ($vec-ptr:(double *qMatMut))[0 * $(int nTs) + j] = NV_Ith_S(y,j);
                         }

                         /* Explicitly set the method */
                         if ($(int method) >= MIN_DIRK_NUM) {
                           /* Implicit */
                           flag = ARKStepSetTableNum(arkode_mem, $(int method), -1);
                         } else {
                           /* Explicit */
                           flag = ARKStepSetTableNum(arkode_mem, -1, $(int method));
                         }
                         if (check_flag(&flag, "ARKStepSetTableNum", 1, report_error)) return 1;

                         /* Main time-stepping loop: calls ARKStepEvolve to perform the integration */
                         /* Stops when the final time has been reached                       */
                         for (i = 1; i < $(int nTs); i++) {

                           flag = ARKStepEvolve(arkode_mem, ($vec-ptr:(double *ts))[i], y, &t, ARK_NORMAL); /* call integrator */
                           if (check_flag(&flag, "ARKode", 1, report_error)) return 1;

                           /* Store the results for Haskell */
                           for (j = 0; j < NEQ; j++) {
                             ($vec-ptr:(double *qMatMut))[i * NEQ + j] = NV_Ith_S(y,j);
                           }
                         }

                         /* Get some final statistics on how the solve progressed */

                         flag = ARKStepGetNumSteps(arkode_mem, &nst);
                         check_flag(&flag, "ARKStepGetNumSteps", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[0] = nst;

                         flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
                         check_flag(&flag, "ARKStepGetNumStepAttempts", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[1] = nst_a;

                         flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
                         check_flag(&flag, "ARKStepGetNumRhsEvals", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[2] = nfe;
                         ($vec-ptr:(sunindextype *diagMut))[3] = nfi;

                         flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
                         check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[4] = nsetups;

                         flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
                         check_flag(&flag, "ARKStepGetNumErrTestFails", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[5] = netf;

                         flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
                         check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[6] = nni;

                         flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
                         check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[7] = ncfn;

                         flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
                         check_flag(&flag, "ARKStepGetNumJacEvals", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[8] = ncfn;

                         flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
                         check_flag(&flag, "ARKStepGetNumRhsEvals", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[9] = ncfn;

                         /* Clean up and return */
                         N_VDestroy(y);            /* Free y vector          */
                         N_VDestroy(tv);           /* Free tv vector         */
                         ARKStepFree(&arkode_mem);  /* Free integrator memory */
                         SUNLinSolFree(LS);        /* Free linear solver     */
                         SUNMatDestroy(A);         /* Free A matrix          */

                         return flag;
                       } |]
  preD <- V.freeze diagMut
  let d = SundialsDiagnostics (fromIntegral $ preD V.!0)
                              (fromIntegral $ preD V.!1)
                              (fromIntegral $ preD V.!2)
                              (fromIntegral $ preD V.!3)
                              (fromIntegral $ preD V.!4)
                              (fromIntegral $ preD V.!5)
                              (fromIntegral $ preD V.!6)
                              (fromIntegral $ preD V.!7)
                              (fromIntegral $ preD V.!8)
                              (fromIntegral $ preD V.!9)
                              False
  m <- V.freeze qMatMut
  if res == 0
    then do
      return $ Right (m, d)
    else do
      return $ Left  (m, res)
