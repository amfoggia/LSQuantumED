#ifndef __PRINT_ARGS_H
#define __TOOLS_H

#include "meson_config.h"
#include <petscsys.h>
#include <mpi.h>

/* --------------------------------------------------------------------------- */
// -------------------------- Printing arguments ----------------------------- //
/* --------------------------------------------------------------------------- */

#ifdef DISORDER
#define print_args()							\
  ierr = PetscPrintf(PETSC_COMM_WORLD, " ------ Parameters ------\n"); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "J1", J1, "D1", Delta1); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "J2", J2, "D2", Delta2); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "min dis", min_dis, "max dis", max_dis); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %d\n%9s:   %d\n", "reprod?", reprod, "repet", rep); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, " ------------------------\n\n"); CHKERRQ(ierr);
#else
#define print_args()							\
  ierr = PetscPrintf(PETSC_COMM_WORLD, " ------ Parameters ------\n"); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "J1", J1, "D1", Delta1); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "J2", J2, "D2", Delta2); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %d\n", "disorder?", false); CHKERRQ(ierr); \
    ierr = PetscPrintf(PETSC_COMM_WORLD, " ------------------------\n\n"); CHKERRQ(ierr);
#endif

#endif
