#ifndef __SOLVER_H
#define __SOLVER_H

#include "meson_config.h"
#include "environment.hpp"
#include <slepceps.h>

/**
 * @namespace Solver
 * @brief Contains solver settings and eigenpairs retrieving functions.
 */

namespace Solver {

  /**
   * @fn Solver::SolverInit(EPS&, Mat, PetscInt, PetscInt, PetscInt)
   * @brief Creates and sets the PETSc EPS object.
   * Sets the type to HERMITIAN, the eigenvalue that it searches to the SMALLEST REAL,
   * and the solver type to KRYLOV-SCHUR. Anyway, allows to overwrite the settings at runtime.
   * @param [in] solver EPS object that contains the solver.
   * @param [in] matrix Hamiltonian matrix.
   * @param [in] nev Number of eigenvalues to compute.
   * @param [in] ncv Maximum dimension of the subspace to be used by the solver.
   * @param [in] mpd Maximum dimension allowed for the projected problem.
   * @return Error value.
   */
  PetscErrorCode SolverInit(EPS& solver,
			    Mat matrix,
			    PetscInt nev,
			    PetscInt ncv = PETSC_DEFAULT,
			    PetscInt mpd = PETSC_DEFAULT);
  
  /**
   * @fn Solver::SolverClean(EPS&)
   * @brief Destroys the PETSc EPS object.
   * @param [in] solver EPS object that contains the solver.
   * @return Error value.
   */
  PetscErrorCode SolverClean(EPS& solver);

  /**
   * @fn Solver::solve(Environment&, EPS&, PetscInt&)
   * @brief Calls PETSc functions to solve the system H x = a x.
   * Gets the number of converged eigenpairs.
   * @param [in] env Environment object.
   * @param [in] solver EPS object that contains the solver.
   * @param [out] nconv Number of converged eigenpairs.
   * @return Error value.
   */
  PetscErrorCode solve(Environment& env, EPS& solver, PetscInt& nconv);

  /**
   * @fn Solver::solution(EPS&, PetscInt, PetscScalar&, Vec&, PetscReal&)
   * @brief Gets the pairs {eigenvalue,eigenvector}.
   * @param [in] solver EPS object that contains the solver.
   * @param [in] index Eigenpair index.
   * @param [in,out] real_value Eigenvalue.
   * @param [in,out] real_vector Eigenvector.
   * @param [in,out] error Error (unless set differently at runtime, relative to the eigenvalue).
   * @return Error value.
   */
  PetscErrorCode solution(EPS& solver,
			  PetscInt index,
			  PetscScalar& real_value,
			  Vec& real_vector,
			  PetscReal& error);
}

#endif
