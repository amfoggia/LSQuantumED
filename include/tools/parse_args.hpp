#ifndef __PARSE_ARGS_H
#define __PARSE_ARGS_H

#include "meson_config.h"
#include <petscsys.h>
#include <mpi.h>

/* --------------------------------------------------------------------------- */
// -------------------------- Parsing arguments ------------------------------ //
/* --------------------------------------------------------------------------- */

// Macro to obtain command-line options for the Hamiltonian
#define parse_args()							\
  PetscInt  _n;								\
  (env.nspins == 0) ? (_n = 6) : (_n = env.nspins);			\
  PetscInt  lx        = 1;						\
  PetscInt  ly        = _n;						\
  PetscReal J1        = 1.0;						\
  PetscReal Delta1    = 1.0;						\
  PetscReal J2        = 0.0;						\
  PetscReal Delta2    = 0.0;						\
  PetscInt  rep       = 1;						\
  PetscReal min_dis   = -0.3;						\
  PetscReal max_dis   = 0.3;						\
  PetscBool reprod    = PETSC_TRUE;					\
  PetscBool check_dis = PETSC_FALSE;					\
  PetscBool check_lattype = PETSC_FALSE;				\
  lattice_type ltype  = lattice_type::chain1D;				\
  const char * const LatTypes[] = {"chain1D", "square2D", "honeycomb2D", "lattice_type", "", 0}; \
									\
									\
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the spin system", ""); CHKERRQ(ierr); \
  {									\
    ierr = PetscOptionsInt("-n", "Number of spins in the system", "environment.hpp", _n, &_n, NULL); CHKERRQ(ierr); \
    if (_n%2 != 0)							\
      throw std::invalid_argument("The total number of spins has to be EVEN.");	\
    if (_n <= 0)							\
      throw std::invalid_argument("The total number of spins has to be POSITIVE."); \
    env.nspins = _n;							\
									\
    ierr = PetscOptionsEnum("-ltype", "Type of lattice", "lattice.hpp", LatTypes, (PetscEnum)ltype, (PetscEnum*)&ltype, &check_lattype); CHKERRQ(ierr); \
    if (ltype == lattice_type::chain1D) {				\
      ierr = PetscOptionsReject(NULL, NULL, "-lx", "For chain1D lattice, one should not specify the dimension in each direction"); CHKERRQ(ierr); \
      ierr = PetscOptionsReject(NULL, NULL, "-ly", "For chain1D lattice, one should not specify the dimension in each direction"); CHKERRQ(ierr); \
    }									\
    else {								\
      ierr = PetscOptionsInt("-lx", "Lattice dimention in x direction", "lattice.hpp", lx, &lx, NULL); CHKERRQ(ierr); \
      ierr = PetscOptionsInt("-ly", "Lattice dimention in y direction", "lattice.hpp", ly, &ly, NULL); CHKERRQ(ierr); \
    }									\
  }									\
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);				\
    									\
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Options for the Hamiltonian", "hamiltonian.hpp"); CHKERRQ(ierr); \
  {									\
    ierr = PetscOptionsReal("-j1", "Nearest neighbours interaction", "Hamiltonian::build_off_diag", J1, &J1, NULL); CHKERRQ(ierr); \
    ierr = PetscOptionsReal("-d1", "Nearest neighbours anisotropy", "Hamiltonian::build_diag", Delta1, &Delta1, NULL); CHKERRQ(ierr); \
    ierr = PetscOptionsReal("-j2", "Next-nearest neighbours interaction", "Hamiltonian::build_off_diag", J2, &J2, NULL); CHKERRQ(ierr); \
    ierr = PetscOptionsReal("-d2", "Next-nearest neighbours anisotropy", "Hamiltonian::build_diag", Delta2, &Delta2, NULL); CHKERRQ(ierr); \
									\
    ierr = PetscOptionsName("-disorder", "Specify if you want to run with disorder average", "Hamiltonian::build_disorder", &check_dis); CHKERRQ(ierr); \
    if (check_dis) {							\
      env.disorder_flg = PETSC_TRUE;					\
      ierr = PetscOptionsReal("-min_dis", "Minimum value for random disorder variable", "Tools::RandomDisorder", min_dis, &min_dis, NULL); CHKERRQ(ierr); \
      ierr = PetscOptionsReal("-max_dis", "Maximum value for random disorder variable", "Tools::RandomDisorder", max_dis, &max_dis, NULL); CHKERRQ(ierr); \
      ierr = PetscOptionsInt("-rep", "Number of disorder repetitions", "Tools::RandomDisorder", rep, &rep, NULL); CHKERRQ(ierr); \
      ierr = PetscOptionsBool("-reprod", "Specify if you want that the disorder random values are reproducible or not", "Tools::RandomDisorder", reprod, &reprod, NULL); CHKERRQ(ierr); \
    }									\
  }									\
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);				\
  check_lattype = PETSC_TRUE;						\
  if (env.disorder_flg == PETSC_FALSE) {				\
    max_dis = 0.0;							\
    min_dis = 0.0;							\
    reprod = PETSC_TRUE;						\
    rep = 1;								\
  }									\

#endif
