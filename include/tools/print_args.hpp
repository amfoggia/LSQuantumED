#ifndef __PRINT_ARGS_H
#define __TOOLS_H

#include "meson_config.h"
#include <petscsys.h>
#include <mpi.h>

/* --------------------------------------------------------------------------- */
// -------------------------- Printing arguments ----------------------------- //
/* --------------------------------------------------------------------------- */

#define print_args()							\
  ierr = PetscPrintf(PETSC_COMM_WORLD, " ------ Parameters ------\n"); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %d\n", "spins", env.nspins); CHKERRQ(ierr); \
  if (check_lattype == true) {						\
    switch(ltype) {							\
    case lattice_type::chain1D: { ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %7s\n", "lattice", "chain1D"); CHKERRQ(ierr); break; } \
    case lattice_type::square2D: {					\
      ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %8s\n", "lattice", "square2D"); CHKERRQ(ierr); \
      ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %d\n%9s:   %d\n", "Lx", lx, "Ly", ly); CHKERRQ(ierr); \
      break;								\
    }									\
    case lattice_type::honeycomb2D: { ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %10s\n", "lattice", "honeycomb2D"); CHKERRQ(ierr); break; } \
    }									\
  }									\
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "J1", J1, "D1", Delta1); CHKERRQ(ierr); \
  ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "J2", J2, "D2", Delta2); CHKERRQ(ierr); \
  if (check_dis == true) {						\
    ierr = PetscPrintf(PETSC_COMM_WORLD, " ------------------------\n"); CHKERRQ(ierr); \
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %4s\n", "disorder", "TRUE"); CHKERRQ(ierr); \
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %.5f\n%9s:   %.5f\n", "min dis", min_dis, "max dis", max_dis); CHKERRQ(ierr); \
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %d\n%9s:   %d\n", "reprod?", reprod, "rep", rep); CHKERRQ(ierr); \
  }									\
  else {								\
    ierr = PetscPrintf(PETSC_COMM_WORLD, " ------------------------\n"); CHKERRQ(ierr); \
    ierr = PetscPrintf(PETSC_COMM_WORLD, "%9s:   %4s\n", "disorder", "FALSE"); CHKERRQ(ierr); \
  }									\
  ierr = PetscPrintf(PETSC_COMM_WORLD, " ------------------------\n\n"); CHKERRQ(ierr);
#endif
