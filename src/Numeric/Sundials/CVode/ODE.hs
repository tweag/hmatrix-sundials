{-# OPTIONS_GHC -Wall -Wno-partial-type-signatures #-}

{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE MultiWayIf #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE NamedFieldPuns #-}

-- | Solution of ordinary differential equation (ODE) initial value problems.
--
-- <https://computation.llnl.gov/projects/sundials/sundials-software>
module Numeric.Sundials.CVode.ODE
  ( odeSolveWithEvents
  , ODEMethod(..)
  , SolverResult(..)
  , cvOdeC
  ) where

import qualified Language.C.Inline as C
import qualified Language.C.Inline.Unsafe as CU

import           Data.Monoid ((<>))
import           Data.Maybe
import           Data.List (genericLength)
import           Control.Applicative

import           Foreign.C.Types (CDouble, CInt)
import           Foreign.Ptr
import           Foreign.Storable (peek, poke)

import qualified Data.Vector.Storable as V
import qualified Data.Vector.Storable.Mutable as VM
import qualified Data.Vector as VB -- B for Boxed

import           Data.Coerce (coerce)

import           Numeric.LinearAlgebra.Devel (createVector)

import           Numeric.LinearAlgebra.HMatrix (Vector, Matrix,
                                                reshape,
                                                subVector, toColumns, fromColumns, asColumn)

import           Numeric.Sundials.Foreign (cV_ADAMS, cV_BDF,
                                          vectorToC, cV_SUCCESS,
                                          SunVector(..), SunIndexType)
import qualified Numeric.Sundials.Foreign as T
import           Numeric.Sundials.Types
import Numeric.Sundials.Common

import Control.Monad.IO.Class
import Katip
import GHC.Prim


C.context (C.baseCtx <> C.vecCtx <> C.funCtx <> sunCtx)

C.include "<stdlib.h>"
C.include "<stdio.h>"
C.include "<string.h>"
C.include "<math.h>"
C.include "<arkode/arkode.h>"
C.include "<cvode/cvode.h>"               -- prototypes for CVODE fcts., consts.
C.include "<nvector/nvector_serial.h>"    -- serial N_Vector types, fcts., macros
C.include "<sunmatrix/sunmatrix_dense.h>" -- access to dense SUNMatrix
C.include "<sunlinsol/sunlinsol_dense.h>" -- access to dense SUNLinearSolver
C.include "<cvode/cvode_direct.h>"        -- access to CVDls interface
C.include "<sundials/sundials_types.h>"   -- definition of type realtype
C.include "<sundials/sundials_math.h>"
C.include "../../../helpers.h"
C.include "Numeric/Sundials/Foreign_hsc.h"

-- | Stepping functions
data ODEMethod = ADAMS
               | BDF
  deriving (Eq, Ord, Show, Read)

instance Method ODEMethod where
  methodToInt ADAMS = cV_ADAMS
  methodToInt BDF   = cV_BDF

cvOdeC :: CConsts -> CVars (V.MVector RealWorld) -> ReportErrorFn -> IO CInt
cvOdeC CConsts{..} CVars{..} report_error =
  [C.block| int {
  /* general problem variables */

  int flag;                  /* reusable error-checking flag                 */

  int i, j;                  /* reusable loop indices                        */
  N_Vector y = NULL;         /* empty vector for storing solution            */
  N_Vector tv = NULL;        /* empty vector for storing absolute tolerances */

  SUNMatrix A = NULL;        /* empty matrix for linear solver               */
  SUNLinearSolver LS = NULL; /* empty linear solver object                   */
  void *cvode_mem = NULL;    /* empty CVODE memory structure                 */
  realtype t;
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;

  realtype tout;

  /* input_ind tracks the current index into the c_sol_time array */
  int input_ind = 1;
  /* output_ind tracks the current row into the c_output_mat matrix.
     If differs from input_ind because of the extra rows corresponding to events. */
  int output_ind = 1;
  /* We need to update c_n_rows every time we update output_ind because
     of the possibility of early return (in which case we still need to assemble
     the partial results matrix). We could even work with c_n_rows only and ditch
     output_ind, but the inline-c expression is quite verbose, and output_ind is
     more convenient to use in index calculations.
  */
  ($vec-ptr:(int *c_n_rows))[0] = output_ind;
  /* event_ind tracks the current event number */
  int event_ind = 0;

  /* general problem parameters */

  realtype T0 = RCONST(($vec-ptr:(double *c_sol_time))[0]); /* initial time              */
  sunindextype c_dim = $(sunindextype c_dim);           /* number of dependent vars. */

  /* Initialize data structures */

  ARKErrHandlerFn report_error = $fun:(void (*report_error)(int,const char*, const char*, char*, void*));

  /* Initialize odeMaxEventsReached to False */
  ($vec-ptr:(sunindextype *c_diagnostics))[10] = 0;

  y = N_VNew_Serial(c_dim); /* Create serial vector for solution */
  if (check_flag((void *)y, "N_VNew_Serial", 0, report_error)) return 1;
  /* Specify initial condition */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(y,i) = ($vec-ptr:(double *c_init_cond))[i];
  };

  // NB: Uses the Newton solver by default
  cvode_mem = CVodeCreate($(int c_method));
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0, report_error)) return(1);

  flag = CVodeSetErrHandlerFn(cvode_mem, report_error, NULL);
  if (check_flag(&flag, "CVodeSetErrHandlerFn", 1, report_error)) return 1;

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, $(int (* c_rhs) (double t, SunVector y[], SunVector dydt[], UserData* params)), T0, y);
  if (check_flag(&flag, "CVodeInit", 1, report_error)) return(1);
  flag = CVodeSetUserData(cvode_mem, $(UserData* c_rhs_userdata));
  if (check_flag(&flag, "CVodeSetUserData", 1, report_error)) return(1);

  tv = N_VNew_Serial(c_dim); /* Create serial vector for absolute tolerances */
  if (check_flag((void *)tv, "N_VNew_Serial", 0, report_error)) return 1;
  /* Specify tolerances */
  for (i = 0; i < c_dim; i++) {
    NV_Ith_S(tv,i) = ($vec-ptr:(double *c_atol))[i];
  };

  flag = CVodeSetMinStep(cvode_mem, $(double c_minstep));
  if (check_flag(&flag, "CVodeSetMinStep", 1, report_error)) return 1;
  flag = CVodeSetMaxNumSteps(cvode_mem, $(sunindextype c_max_n_steps));
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1, report_error)) return 1;
  flag = CVodeSetMaxErrTestFails(cvode_mem, $(int c_max_err_test_fails));
  if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1, report_error)) return 1;

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, $(double c_rtol), tv);
  if (check_flag(&flag, "CVodeSVtolerances", 1, report_error)) return(1);

  /* Call CVodeRootInit to specify the root function c_event_fn with c_n_event_specs components */
  flag = CVodeRootInit(cvode_mem, $(int c_n_event_specs), $fun:(int (* c_event_fn) (double t, SunVector y[], double gout[], void * params)));

  if (check_flag(&flag, "CVodeRootInit", 1, report_error)) return(1);

  /* Initialize dense matrix data structure and solver */
  A = SUNDenseMatrix(c_dim, c_dim);
  if (check_flag((void *)A, "SUNDenseMatrix", 0, report_error)) return 1;
  LS = SUNDenseLinearSolver(y, A);
  if (check_flag((void *)LS, "SUNDenseLinearSolver", 0, report_error)) return 1;

  /* Attach matrix and linear solver */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if (check_flag(&flag, "CVDlsSetLinearSolver", 1, report_error)) return 1;

  /* Set the initial step size if there is one */
  if ($(int c_init_step_size_set)) {
    /* FIXME: We could check if the initial step size is 0 */
    /* or even NaN and then throw an error                 */
    flag = CVodeSetInitStep(cvode_mem, $(double c_init_step_size));
    if (check_flag(&flag, "CVodeSetInitStep", 1, report_error)) return 1;
  }

  /* Set the Jacobian if there is one */
  if ($(int c_jac_set)) {
    flag = CVDlsSetJacFn(cvode_mem, $fun:(int (* c_jac) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
    if (check_flag(&flag, "CVDlsSetJacFn", 1, report_error)) return 1;
  }

  /* Store initial conditions */
  ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + 0] = ($vec-ptr:(double *c_sol_time))[0];
  for (j = 0; j < c_dim; j++) {
    ($vec-ptr:(double *c_output_mat))[0 * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
  }

  while (1) {
    flag = CVode(cvode_mem, ($vec-ptr:(double *c_sol_time))[input_ind], y, &t, CV_NORMAL); /* call integrator */
    if (check_flag(&flag, "CVode", 1, report_error)) {
      N_Vector ele = N_VNew_Serial(c_dim);
      N_Vector weights = N_VNew_Serial(c_dim);
      flag = CVodeGetEstLocalErrors(cvode_mem, ele);
      // CV_SUCCESS is defined is 0, so we OR the flags
      flag = flag || CVodeGetErrWeights(cvode_mem, weights);
      if (flag == CV_SUCCESS) {
        double *arr_ptr = N_VGetArrayPointer(ele);
        memcpy(($vec-ptr:(double *c_local_error)), arr_ptr, c_dim * sizeof(double));

        arr_ptr = N_VGetArrayPointer(weights);
        memcpy(($vec-ptr:(double *c_var_weight)), arr_ptr, c_dim * sizeof(double));

        ($vec-ptr:(int *c_local_error_set))[0] = 1;
      }
      N_VDestroy(ele);
      N_VDestroy(weights);
      return 1;
    }

    /* Store the results for Haskell */
    ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
    for (j = 0; j < c_dim; j++) {
      ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
    }
    output_ind++;
    ($vec-ptr:(int *c_n_rows))[0] = output_ind;

    if (flag == CV_ROOT_RETURN) {
      if (event_ind >= $(int c_max_events)) {
        /* We reached the maximum number of events.
           Either the maximum number of events is set to 0,
           or there's a bug in our code below. In any case return an error.
        */
        return 1;
      }

      /* Are we interested in this event?
         If not, continue without any observable side-effects.
      */
      int good_event = 0;
      int stop_solver = 0;
      flag = CVodeGetRootInfo(cvode_mem, ($vec-ptr:(int *c_root_info)));
      if (check_flag(&flag, "CVodeGetRootInfo", 1, report_error)) return 1;
      for (i = 0; i < $(int c_n_event_specs); i++) {
        int ev = ($vec-ptr:(int *c_root_info))[i];
        int req_dir = ($vec-ptr:(const int *c_requested_event_direction))[i];
        if (ev != 0 && ev * req_dir >= 0) {
          good_event = 1;

          ($vec-ptr:(int *c_actual_event_direction))[event_ind] = ev;
          ($vec-ptr:(int *c_event_index))[event_ind] = i;
          ($vec-ptr:(double *c_event_time))[event_ind] = t;
          event_ind++;
          stop_solver = ($vec-ptr:(int *c_event_stops_solver))[i];

          /* Update the state with the supplied function */
          $fun:(int (* c_apply_event) (int, double, SunVector y[], SunVector z[]))(i, t, y, y);
        }
      }

      if (good_event) {
        ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + 0] = t;
        for (j = 0; j < c_dim; j++) {
          ($vec-ptr:(double *c_output_mat))[output_ind * (c_dim + 1) + (j + 1)] = NV_Ith_S(y,j);
        }
        output_ind++;
        ($vec-ptr:(int *c_n_rows))[0] = output_ind;

        if (stop_solver) {
          break;
        }
        if (event_ind >= $(int c_max_events)) {
          /* We collected the requested number of events. Stop the solver. */
          ($vec-ptr:(sunindextype *c_diagnostics))[10] = 1;
          break;
        }

        flag = CVodeReInit(cvode_mem, t, y);
        if (check_flag(&flag, "CVodeReInit", 1, report_error)) return(1);
      } else {
        /* Since this is not a wanted event, it shouldn't get a row */
        output_ind--;
        ($vec-ptr:(int *c_n_rows))[0] = output_ind;
      }
    }
    else {
      if (++input_ind >= $(int c_n_sol_times))
        break;
    }
  }

  /* The number of actual roots we found */
  ($vec-ptr:(int *c_n_events))[0] = event_ind;

  /* Get some final statistics on how the solve progressed */
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[0] = nst;

  /* FIXME */
  ($vec-ptr:(sunindextype *c_diagnostics))[1] = 0;

  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[2] = nfe;
  /* FIXME */
  ($vec-ptr:(sunindextype *c_diagnostics))[3] = 0;

  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[4] = nsetups;

  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[5] = netf;

  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[6] = nni;

  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[7] = ncfn;

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[8] = ncfn;

  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1, report_error);
  ($vec-ptr:(sunindextype *c_diagnostics))[9] = ncfn;

  /* Clean up and return */

  N_VDestroy(y);          /* Free y vector          */
  N_VDestroy(tv);         /* Free tv vector         */
  CVodeFree(&cvode_mem);  /* Free integrator memory */
  SUNLinSolFree(LS);      /* Free linear solver     */
  SUNMatDestroy(A);       /* Free A matrix          */

  return CV_SUCCESS;
 } |]

solveOdeC
  :: (Katip m, MonadIO m)
  => CInt
  -> SunIndexType
  -> CDouble
  -> CInt
  -> Maybe CDouble
  -> (Maybe (CDouble -> V.Vector CDouble -> T.SunMatrix))
  -> Tolerances
  -> OdeRhs -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> V.Vector CDouble -- ^ Initial conditions
  -> CInt -- ^ Number of event equations
  -> (CDouble -> V.Vector CDouble -> V.Vector CDouble) -- ^ The event equations themselves
  -> VB.Vector CrossingDirection -- ^ The required crossing direction for each event
  -> V.Vector CInt -- ^ Whether an event should stop the solver
  -> CInt -- ^ Maximum number of events
  -> (Int -> CDouble -> V.Vector CDouble -> V.Vector CDouble)
      -- ^ Function to reset/update the state when an event occurs. The
      -- 'Int' argument is the 0-based number of the event that has
      -- occurred. If multiple events have occurred at the same time, they
      -- are handled in the increasing order of the event index. The other
      -- arguments are the time and the point in the state space. Return
      -- the updated point in the state space.
  -> V.Vector CDouble -- ^ Desired solution times
  -> m SolverResult
solveOdeC maxErrTestFails maxNumSteps_ minStep_ method initStepSize
          jacH (Tolerances rTol aTols0) rhs f0 nr event_fn directions event_stops_solver max_events apply_event ts
  | V.null f0 = -- 0-dimensional (empty) system
    return $ SolverSuccess [] (asColumn (coerce ts)) emptyDiagnostics
  | otherwise = do

  report_error <- logWithKatip

  liftIO $ do -- the rest is in the IO monad

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
      nEq :: SunIndexType
      nEq = fromIntegral dim
      nTs :: CInt
      nTs = fromIntegral $ V.length ts
      aTols :: V.Vector CDouble
      aTols = either (V.replicate dim) id aTols0

  output_mat_mut :: V.MVector _ CDouble <- V.thaw =<< createVector ((1 + fromIntegral dim) * (fromIntegral (2 * max_events) + fromIntegral nTs))
  -- diagMut is a mutable vector which we write diagnostic data while
  -- solving. Its size corresponds to the number of fields in
  -- SundialsDiagnostics.
  diagMut :: V.MVector _ SunIndexType <- V.thaw =<< createVector 11 -- FIXME
  (rhs_funptr :: FunPtr OdeRhsCType, userdata :: Ptr UserData) <-
    case rhs of
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

  let nrPre = fromIntegral nr
  gResults :: V.Vector CInt <- createVector nrPre
  -- FIXME: Do we need to do this here? Maybe as it will get GC'd and
  -- we'd have to do a malloc in C otherwise :(
  gResMut <- V.thaw gResults
  event_index_mut :: V.MVector _ CInt <- V.thaw =<< createVector (fromIntegral max_events)
  event_time_mut :: V.MVector _ CDouble <- V.thaw =<< createVector (fromIntegral max_events)
  -- Total number of events. This is *not* directly re
  n_events_mut :: V.MVector _ CInt <- V.thaw =<< createVector 1
  -- Total number of rows in the output_mat_mut matrix. It *cannot* be
  -- inferred from n_events_mut because when an event occurs k times, it
  -- contributes k to n_events_mut but only 2 to n_rows_mut.
  n_rows_mut :: V.MVector _ CInt <- V.thaw =<< createVector 1
  actual_event_direction_mut :: V.MVector _ CInt <- V.thaw =<< createVector (fromIntegral max_events)

  -- The vector containing local error estimates
  local_errors_mut :: V.MVector _ CDouble <- V.thaw =<< createVector dim
  -- The vector containing variable weights
  var_weights_mut :: V.MVector _ CDouble <- V.thaw =<< createVector dim
  -- The flag indicating whether local_errors_mut is filled with meaningful
  -- values
  local_errors_set :: V.MVector _ CInt <- V.thaw =<< createVector 1
  VM.write local_errors_set 0 0

  let event_fn_c :: CDouble -> Ptr T.SunVector -> Ptr CDouble -> Ptr () -> IO CInt
      event_fn_c x y f _ptr = do
        vals <- event_fn x <$> (sunVecVals <$> peek y)
        -- FIXME: We should be able to use poke somehow
        vectorToC vals nrPre f
        return 0

      apply_event_c :: CInt -> CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> IO CInt
      apply_event_c event_index t y y' = do
        sv <- peek y
        poke y' $ SunVector
          { sunVecN = sunVecN sv
          , sunVecVals = apply_event (fromIntegral event_index) t (sunVecVals sv)
          }
        return 0

      requested_event_directions :: V.Vector CInt
      requested_event_directions = V.convert $ VB.map directionToInt directions

      isJac :: CInt
      isJac = fromIntegral $ fromEnum $ isJust jacH
      jacIO :: CDouble -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunMatrix ->
               Ptr () -> Ptr T.SunVector -> Ptr T.SunVector -> Ptr T.SunVector ->
               IO CInt
      jacIO t y _fy jacS _ptr _tmp1 _tmp2 _tmp3 = do
        case jacH of
          Nothing   -> error "Numeric.Sundials.CVode.ODE: Jacobian not defined"
          Just jacI -> do j <- jacI t <$> (sunVecVals <$> peek y)
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
                         void *cvode_mem = NULL;    /* empty CVODE memory structure                 */
                         realtype t;
                         long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;

                         realtype tout;

                         /* input_ind tracks the current index into the ts array */
                         int input_ind = 1;
                         /* output_ind tracks the current row into the output_mat_mut matrix.
                            If differs from input_ind because of the extra rows corresponding to events. */
                         int output_ind = 1;
                         /* We need to update n_rows_mut every time we update output_ind because
                            of the possibility of early return (in which case we still need to assemble
                            the partial results matrix). We could even work with n_rows_mut only and ditch
                            output_ind, but the inline-c expression is quite verbose, and output_ind is
                            more convenient to use in index calculations.
                         */
                         ($vec-ptr:(int *n_rows_mut))[0] = output_ind;
                         /* event_ind tracks the current event number */
                         int event_ind = 0;

                         /* general problem parameters */

                         realtype T0 = RCONST(($vec-ptr:(double *ts))[0]); /* initial time              */
                         sunindextype NEQ = $(sunindextype nEq);           /* number of dependent vars. */

                         /* Initialize data structures */

                         ARKErrHandlerFn report_error = $fun:(void (*report_error)(int,const char*, const char*, char*, void*));

                         /* Initialize odeMaxEventsReached to False */
                         ($vec-ptr:(sunindextype *diagMut))[10] = 0;

                         y = N_VNew_Serial(NEQ); /* Create serial vector for solution */
                         if (check_flag((void *)y, "N_VNew_Serial", 0, report_error)) return 1;
                         /* Specify initial condition */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(y,i) = ($vec-ptr:(double *f0))[i];
                         };

                         // NB: Uses the Newton solver by default
                         cvode_mem = CVodeCreate($(int method));
                         if (check_flag((void *)cvode_mem, "CVodeCreate", 0, report_error)) return(1);

                         flag = CVodeSetErrHandlerFn(cvode_mem, report_error, NULL);
                         if (check_flag(&flag, "CVodeSetErrHandlerFn", 1, report_error)) return 1;

                         /* Call CVodeInit to initialize the integrator memory and specify the
                          * user's right hand side function in y'=f(t,y), the inital time T0, and
                          * the initial dependent variable vector y. */
                         flag = CVodeInit(cvode_mem, $(int (* rhs_funptr) (double t, SunVector y[], SunVector dydt[], UserData* params)), T0, y);
                         if (check_flag(&flag, "CVodeInit", 1, report_error)) return(1);
                         flag = CVodeSetUserData(cvode_mem, $(UserData* userdata));
                         if (check_flag(&flag, "CVodeSetUserData", 1, report_error)) return(1);

                         tv = N_VNew_Serial(NEQ); /* Create serial vector for absolute tolerances */
                         if (check_flag((void *)tv, "N_VNew_Serial", 0, report_error)) return 1;
                         /* Specify tolerances */
                         for (i = 0; i < NEQ; i++) {
                           NV_Ith_S(tv,i) = ($vec-ptr:(double *aTols))[i];
                         };

                         flag = CVodeSetMinStep(cvode_mem, $(double minStep_));
                         if (check_flag(&flag, "CVodeSetMinStep", 1, report_error)) return 1;
                         flag = CVodeSetMaxNumSteps(cvode_mem, $(sunindextype maxNumSteps_));
                         if (check_flag(&flag, "CVodeSetMaxNumSteps", 1, report_error)) return 1;
                         flag = CVodeSetMaxErrTestFails(cvode_mem, $(int maxErrTestFails));
                         if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1, report_error)) return 1;

                         /* Call CVodeSVtolerances to specify the scalar relative tolerance
                          * and vector absolute tolerances */
                         flag = CVodeSVtolerances(cvode_mem, $(double rTol), tv);
                         if (check_flag(&flag, "CVodeSVtolerances", 1, report_error)) return(1);

                         /* Call CVodeRootInit to specify the root function event_fn_c with nr components */
                         flag = CVodeRootInit(cvode_mem, $(int nr), $fun:(int (* event_fn_c) (double t, SunVector y[], double gout[], void * params)));

                         if (check_flag(&flag, "CVodeRootInit", 1, report_error)) return(1);

                         /* Initialize dense matrix data structure and solver */
                         A = SUNDenseMatrix(NEQ, NEQ);
                         if (check_flag((void *)A, "SUNDenseMatrix", 0, report_error)) return 1;
                         LS = SUNDenseLinearSolver(y, A);
                         if (check_flag((void *)LS, "SUNDenseLinearSolver", 0, report_error)) return 1;

                         /* Attach matrix and linear solver */
                         flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
                         if (check_flag(&flag, "CVDlsSetLinearSolver", 1, report_error)) return 1;

                         /* Set the initial step size if there is one */
                         if ($(int isInitStepSize)) {
                           /* FIXME: We could check if the initial step size is 0 */
                           /* or even NaN and then throw an error                 */
                           flag = CVodeSetInitStep(cvode_mem, $(double ss));
                           if (check_flag(&flag, "CVodeSetInitStep", 1, report_error)) return 1;
                         }

                         /* Set the Jacobian if there is one */
                         if ($(int isJac)) {
                           flag = CVDlsSetJacFn(cvode_mem, $fun:(int (* jacIO) (double t, SunVector y[], SunVector fy[], SunMatrix Jac[], void * params, SunVector tmp1[], SunVector tmp2[], SunVector tmp3[])));
                           if (check_flag(&flag, "CVDlsSetJacFn", 1, report_error)) return 1;
                         }

                         /* Store initial conditions */
                         ($vec-ptr:(double *output_mat_mut))[0 * (NEQ + 1) + 0] = ($vec-ptr:(double *ts))[0];
                         for (j = 0; j < NEQ; j++) {
                           ($vec-ptr:(double *output_mat_mut))[0 * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
                         }

                         while (1) {
                           flag = CVode(cvode_mem, ($vec-ptr:(double *ts))[input_ind], y, &t, CV_NORMAL); /* call integrator */
                           if (check_flag(&flag, "CVode", 1, report_error)) {
                             N_Vector ele = N_VNew_Serial(NEQ);
                             N_Vector weights = N_VNew_Serial(NEQ);
                             flag = CVodeGetEstLocalErrors(cvode_mem, ele);
                             // CV_SUCCESS is defined is 0, so we OR the flags
                             flag = flag || CVodeGetErrWeights(cvode_mem, weights);
                             if (flag == CV_SUCCESS) {
                               double *arr_ptr = N_VGetArrayPointer(ele);
                               memcpy(($vec-ptr:(double *local_errors_mut)), arr_ptr, NEQ * sizeof(double));

                               arr_ptr = N_VGetArrayPointer(weights);
                               memcpy(($vec-ptr:(double *var_weights_mut)), arr_ptr, NEQ * sizeof(double));

                               ($vec-ptr:(int *local_errors_set))[0] = 1;
                             }
                             N_VDestroy(ele);
                             N_VDestroy(weights);
                             return 1;
                           }

                           /* Store the results for Haskell */
                           ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + 0] = t;
                           for (j = 0; j < NEQ; j++) {
                             ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
                           }
                           output_ind++;
                           ($vec-ptr:(int *n_rows_mut))[0] = output_ind;

                           if (flag == CV_ROOT_RETURN) {
                             if (event_ind >= $(int max_events)) {
                               /* We reached the maximum number of events.
                                  Either the maximum number of events is set to 0,
                                  or there's a bug in our code below. In any case return an error.
                               */
                               return 1;
                             }

                             /* Are we interested in this event?
                                If not, continue without any observable side-effects.
                             */
                             int good_event = 0;
                             int stop_solver = 0;
                             flag = CVodeGetRootInfo(cvode_mem, ($vec-ptr:(int *gResMut)));
                             if (check_flag(&flag, "CVodeGetRootInfo", 1, report_error)) return 1;
                             for (i = 0; i < $(int nr); i++) {
                               int ev = ($vec-ptr:(int *gResMut))[i];
                               int req_dir = ($vec-ptr:(const int *requested_event_directions))[i];
                               if (ev != 0 && ev * req_dir >= 0) {
                                 good_event = 1;

                                 ($vec-ptr:(int *actual_event_direction_mut))[event_ind] = ev;
                                 ($vec-ptr:(int *event_index_mut))[event_ind] = i;
                                 ($vec-ptr:(double *event_time_mut))[event_ind] = t;
                                 event_ind++;
                                 stop_solver = ($vec-ptr:(int *event_stops_solver))[i];

                                 /* Update the state with the supplied function */
                                 $fun:(int (* apply_event_c) (int, double, SunVector y[], SunVector z[]))(i, t, y, y);
                               }
                             }

                             if (good_event) {
                               ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + 0] = t;
                               for (j = 0; j < NEQ; j++) {
                                 ($vec-ptr:(double *output_mat_mut))[output_ind * (NEQ + 1) + (j + 1)] = NV_Ith_S(y,j);
                               }
                               output_ind++;
                               ($vec-ptr:(int *n_rows_mut))[0] = output_ind;

                               if (stop_solver) {
                                 break;
                               }
                               if (event_ind >= $(int max_events)) {
                                 /* We collected the requested number of events. Stop the solver. */
                                 ($vec-ptr:(sunindextype *diagMut))[10] = 1;
                                 break;
                               }

                               flag = CVodeReInit(cvode_mem, t, y);
                               if (check_flag(&flag, "CVodeReInit", 1, report_error)) return(1);
                             } else {
                               /* Since this is not a wanted event, it shouldn't get a row */
                               output_ind--;
                               ($vec-ptr:(int *n_rows_mut))[0] = output_ind;
                             }
                           }
                           else {
                             if (++input_ind >= $(int nTs))
                               break;
                           }
                         }

                         /* The number of actual roots we found */
                         ($vec-ptr:(int *n_events_mut))[0] = event_ind;

                         /* Get some final statistics on how the solve progressed */
                         flag = CVodeGetNumSteps(cvode_mem, &nst);
                         check_flag(&flag, "CVodeGetNumSteps", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[0] = nst;

                         /* FIXME */
                         ($vec-ptr:(sunindextype *diagMut))[1] = 0;

                         flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
                         check_flag(&flag, "CVodeGetNumRhsEvals", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[2] = nfe;
                         /* FIXME */
                         ($vec-ptr:(sunindextype *diagMut))[3] = 0;

                         flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
                         check_flag(&flag, "CVodeGetNumLinSolvSetups", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[4] = nsetups;

                         flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
                         check_flag(&flag, "CVodeGetNumErrTestFails", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[5] = netf;

                         flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
                         check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[6] = nni;

                         flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
                         check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[7] = ncfn;

                         flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
                         check_flag(&flag, "CVDlsGetNumJacEvals", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[8] = ncfn;

                         flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
                         check_flag(&flag, "CVDlsGetNumRhsEvals", 1, report_error);
                         ($vec-ptr:(sunindextype *diagMut))[9] = ncfn;

                         /* Clean up and return */

                         N_VDestroy(y);          /* Free y vector          */
                         N_VDestroy(tv);         /* Free tv vector         */
                         CVodeFree(&cvode_mem);  /* Free integrator memory */
                         SUNLinSolFree(LS);      /* Free linear solver     */
                         SUNMatDestroy(A);       /* Free A matrix          */

                         return CV_SUCCESS;
                       } |]

  -- Free the allocated FunPtr. Ideally this should be done within
  -- a bracket...
  case rhs of
    OdeRhsHaskell {} -> freeHaskellFunPtr rhs_funptr
    OdeRhsC {} -> return () -- we didn't allocate this

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
                              (toEnum . fromIntegral $ preD V.! 10)
  n_rows <- fromIntegral . V.head <$> V.freeze n_rows_mut
  output_mat <- coerce . reshape (dim + 1) . subVector 0 ((dim + 1) * n_rows) <$>
    V.freeze output_mat_mut
  n_events <- fromIntegral . V.head <$> V.freeze n_events_mut
  event_time             :: V.Vector Double
    <- coerce . V.take n_events <$> V.freeze event_time_mut
  event_index            :: V.Vector Int
    <- V.map fromIntegral . V.take n_events <$> V.freeze event_index_mut
  actual_event_direction :: V.Vector CInt
    <- V.take n_events <$> V.freeze actual_event_direction_mut
  (local_errors, var_weights) <- do
    set <- VM.read local_errors_set 0
    if set == 1
      then coerce <$>
        (liftA2 (,)
          (V.freeze local_errors_mut)
          (V.freeze var_weights_mut))
      else mempty

  let
    events :: [EventInfo]
    events = zipWith3 EventInfo
      (V.toList event_time)
      (V.toList event_index)
      (map (fromJust . intToDirection) $ V.toList actual_event_direction)
  return $
    if res == cV_SUCCESS
      then
        SolverSuccess events output_mat d
      else
        SolverError ErrorDiagnostics
          { partialResults = output_mat
          , errorCode = fromIntegral res
          , errorEstimates = local_errors
          , varWeights = var_weights
          }

data SolverResult
  = SolverError !ErrorDiagnostics
  | SolverSuccess
      [EventInfo]
      !(Matrix Double)
      !SundialsDiagnostics
      -- ^ Times at which the event was triggered, information about which root and the
                                                   -- results and diagnostics.
    deriving Show

odeSolveRootVWith'
  :: (MonadIO m, Katip m)
  => ODEOpts ODEMethod
  -> OdeRhs
      -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> Maybe (Double -> Vector Double -> Matrix Double)
      -- ^ The Jacobian (optional)
  -> V.Vector Double                      -- ^ Initial conditions
  -> [EventSpec]                          -- ^ Event specifications
  -> Int                                  -- ^ Maximum number of events
  -> V.Vector Double                      -- ^ Desired solution times
  -> m SolverResult
odeSolveRootVWith' opts rhs mb_jacobian y0 event_specs nRootEvs tt =
  solveOdeC (fromIntegral $ maxFail opts)
                 (fromIntegral $ maxNumSteps opts) (coerce $ minStep opts)
                 (methodToInt . odeMethod $ opts) (coerce $ initStep opts) jacH (stepControl opts)
                 rhs (coerce y0)
                 (genericLength event_specs)
                 event_equations
                 event_directions
                 event_stops_solver
                 (fromIntegral nRootEvs) update_state
                 (coerce tt)
  where
    jacH = fmap (\g t v -> matrixToSunMatrix $ g (coerce t) (coerce v)) $ mb_jacobian
    event_vec :: VB.Vector EventSpec
    event_vec = VB.fromList event_specs
    event_equations :: CDouble -> Vector CDouble -> Vector CDouble
    event_equations t y = V.convert $
      VB.map (\ev -> coerce (eventCondition ev) t y) event_vec
    event_directions :: VB.Vector CrossingDirection
    event_directions = VB.map eventDirection event_vec
    event_stops_solver :: V.Vector CInt
    event_stops_solver =
      V.convert
      . VB.map (fromIntegral . fromEnum . eventStopSolver)
      $ event_vec
    update_state :: Int -> CDouble -> Vector CDouble -> Vector CDouble
    update_state n_event = coerce $ eventUpdate (event_vec VB.! n_event)

odeSolveWithEvents
  :: (MonadIO m, Katip m)
  => ODEOpts ODEMethod
  -> [EventSpec]
    -- ^ Event specifications
  -> Int
    -- ^ Maximum number of events
  -> OdeRhs
    -- ^ The RHS of the system \(\dot{y} = f(t,y)\)
  -> Maybe (Double -> Vector Double -> Matrix Double)
    -- ^ The Jacobian (optional)
  -> V.Vector Double
    -- ^ Initial conditions
  -> V.Vector Double
    -- ^ Desired solution times
  -> m (Either ErrorDiagnostics SundialsSolution)
    -- ^ Either an error code or a solution
odeSolveWithEvents opts event_specs max_events rhs mb_jacobian initial sol_times = do
  result :: SolverResult
    <- odeSolveRootVWith' opts rhs mb_jacobian initial event_specs
        max_events sol_times
  return $ case result of
    SolverError diagn -> Left diagn
    SolverSuccess events mx diagn ->
      Right $ SundialsSolution
          { actualTimeGrid = extractTimeGrid mx
          , solutionMatrix = dropTimeGrid mx
          , eventInfo = events
          , diagnostics = diagn
          }
  where
    -- The time grid is the first column of the result matrix
    extractTimeGrid :: Matrix Double -> Vector Double
    extractTimeGrid = head . toColumns
    dropTimeGrid :: Matrix Double -> Matrix Double
    dropTimeGrid = fromColumns . tail . toColumns
