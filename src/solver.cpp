#include "solver.hpp"

/* --------------------------------------------------------------------------- */
// ----------------------- "Constructors/Destructors" ------------------------ //
/* --------------------------------------------------------------------------- */

PetscErrorCode Solver::SolverInit(EPS& solver,
				  Mat matrix,
				  PetscInt nev,
				  PetscInt ncv,
				  PetscInt mpd) {

  PetscErrorCode ierr = 0;

  ierr = EPSCreate(PETSC_COMM_WORLD, &solver); CHKERRQ(ierr);
  ierr = EPSSetOperators(solver,matrix,NULL); CHKERRQ(ierr);
  ierr = EPSSetProblemType(solver, EPS_HEP); CHKERRQ(ierr);
  ierr = EPSSetType(solver, EPSKRYLOVSCHUR);
  ierr = EPSSetWhichEigenpairs(solver, EPS_SMALLEST_REAL); CHKERRQ(ierr);
  ierr = EPSSetDimensions(solver, nev, ncv, mpd); CHKERRQ(ierr);
  ierr = EPSSetFromOptions(solver); CHKERRQ(ierr);

  return ierr;
}

/* --------------------------------------------------------------------------- */

PetscErrorCode Solver::SolverClean(EPS& solver) {
  PetscErrorCode ierr = 0;
  ierr = EPSDestroy(&solver); CHKERRQ(ierr);
  return ierr;
}

/* --------------------------------------------------------------------------- */
// -------------------------- Solving functions ------------------------------ //
/* --------------------------------------------------------------------------- */

PetscErrorCode Solver::solve(Environment& env, EPS& solver, PetscInt& nconv) {
  
  PetscErrorCode ierr = 0;
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "EPSSolve"};
#endif
    ierr = EPSSolve(solver);
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("EPSSolve");
#endif
  ierr = EPSGetConverged(solver,&nconv);
  return ierr;
}

/* --------------------------------------------------------------------------- */
// -------------------------- Information functions -------------------------- //
/* --------------------------------------------------------------------------- */

PetscErrorCode Solver::solution(EPS& solver,
				PetscInt index,
				PetscScalar& real_value,
				Vec& real_vector,
				PetscReal& error){
  
  PetscErrorCode ierr = 0;
  ierr = EPSGetEigenpair(solver, index, &real_value, NULL, real_vector, NULL); CHKERRQ(ierr);
  ierr = EPSComputeError(solver, index, EPS_ERROR_RELATIVE, &error); CHKERRQ(ierr);
  return ierr;
}
