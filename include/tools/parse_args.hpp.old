#ifndef __PARSE_ARGS_H
#define __PARSE_ARGS_H

#include "meson_config.h"
#include <petscsys.h>
#include <mpi.h>

/* --------------------------------------------------------------------------- */
// -------------------------- Parsing arguments ------------------------------ //
/* --------------------------------------------------------------------------- */

// Macro to obtain command-line options for the Hamiltonian
#ifdef DISORDER
#define parse_args()							\
  PetscReal J1 = 0.0;							\
  PetscReal Delta1 = 0.0;						\
  PetscReal J2 = 0.0;							\
  PetscReal Delta2 = 0.0;						\
  PetscInt rep = 1;							\
  PetscReal min_dis = 0.0, max_dis = 0.0;				\
  PetscBool check_arg_J = PETSC_FALSE,					\
    check_arg_D = PETSC_FALSE,						\
    check_arg_nn = PETSC_FALSE,						\
    check_arg_nnn = PETSC_FALSE;					\
  PetscBool reprod = PETSC_TRUE;					\
  PetscBool check_arg_dis = PETSC_FALSE,				\
    check_param_dis = PETSC_FALSE;					\
									\
  ierr = PetscOptionsHasName(NULL, NULL, "-nn", &check_arg_nn); CHKERRQ(ierr); \
  ierr = PetscOptionsHasName(NULL, NULL, "-nnn", &check_arg_nnn); CHKERRQ(ierr); \
  ierr = PetscOptionsHasName(NULL, NULL, "-disorder", &check_arg_dis); CHKERRQ(ierr); \
									\
  if (!check_arg_nn && !check_arg_nnn && !check_arg_dis) {		\
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide at least one kind of interaction between spins:\n -nn: nearest neighbours\n -nnn: next-nearest neighbours\n -disorder: disorder\n"); CHKERRQ(ierr); \
    return(1);								\
  }									\
  									\
  /*									\
  PetscInt  _latT, Lx, Ly; \
  PetscBool check_arg_lat,					\
    check_arg_lat_x,							\
    check_arg_lat_y;							\
  ierr = PetscOptionsHasName(NULL, NULL, "-lat", &check_arg_lat); CHKERRQ(ierr); \
  if (!check_arg_lat) {							\
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the lattice type with option -lat:\n <0> chain1D, <1> square2D \n"); \
    return 1;								\
  }									\
  if (check_arg_lat) {							\
    ierr = PetscOptionsGetInt(NULL, NULL, "-lat", &_latT, &check_arg_lat); CHKERRQ(ierr); \
    if (_latT < 0 || _latT > 1) {					\
      ierr = PetscPrintf(PETSC_COMM_WORLD, "The allowed options are:\n <0> chain1D, <1> square2D \n"); \
      return 1;								\
    }									\
    if (_latT == 1) {							\
      ierr = PetscOptionsHasName(NULL, NULL, "-Lx", &check_arg_lat_x); CHKERRQ(ierr); \
      ierr = PetscOptionsHasName(NULL, NULL, "-Ly", &check_arg_lat_y); CHKERRQ(ierr); \
      if (check_arg_lat_x && check_arg_lat_y) {				\
	ierr = PetscOptionsGetInt(NULL, NULL, "-Lx", &Lx, &check_arg_lat_x); CHKERRQ(ierr); \
	ierr = PetscOptionsGetInt(NULL, NULL, "-Ly", &Ly, &check_arg_lat_y); CHKERRQ(ierr); \
      }									\
      else {								\
	ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the dimensions of the square lattice with -Lx and -Ly\n"); CHKERRQ(ierr); \
	return 1;							\
      }									\
    }									\
  }									\
  */									\
									\
  if (check_arg_dis) {							\
    ierr = PetscOptionsHasName(NULL, NULL, "-min_dis", &check_param_dis); CHKERRQ(ierr); \
    if (check_param_dis)						\
      ierr = PetscOptionsGetReal(NULL, NULL, "-min_dis", &min_dis, &check_param_dis); CHKERRQ(ierr); \
    ierr = PetscOptionsHasName(NULL, NULL, "-max_dis", &check_param_dis); CHKERRQ(ierr); \
    if (check_param_dis)						\
      ierr = PetscOptionsGetReal(NULL, NULL, "-max_dis", &max_dis, &check_param_dis); CHKERRQ(ierr); \
    ierr = PetscOptionsHasName(NULL, NULL, "-reprod", &check_param_dis); CHKERRQ(ierr); \
    if (check_param_dis)						\
      ierr = PetscOptionsGetBool(NULL, NULL, "-reprod", &reprod, &check_param_dis); CHKERRQ(ierr); \
    ierr = PetscOptionsHasName(NULL, NULL, "-rep", &check_param_dis); CHKERRQ(ierr); \
    if (check_param_dis)						\
      ierr = PetscOptionsGetInt(NULL, NULL, "-rep", &rep, &check_param_dis); CHKERRQ(ierr); \
  }									\
									\
  if (check_arg_nn) {							\
    ierr = PetscOptionsHasName(NULL, NULL, "-J1", &check_arg_J); CHKERRQ(ierr); \
    ierr = PetscOptionsHasName(NULL, NULL, "-D1", &check_arg_D); CHKERRQ(ierr); \
    if (check_arg_J && check_arg_D) {					\
      ierr = PetscOptionsGetReal(NULL, NULL, "-J1", &J1, &check_arg_J); CHKERRQ(ierr); \
      ierr = PetscOptionsGetReal(NULL, NULL, "-D1", &Delta1, &check_arg_D); CHKERRQ(ierr); \
    } else {								\
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the value of the nearest neighbours interaction with option -J1 and anisotropy with option -D1\n"); CHKERRQ(ierr); \
      return(1);							\
    }									\
  }									\
  									\
  if (check_arg_nnn) {							\
    ierr = PetscOptionsHasName(NULL, NULL, "-J2", &check_arg_J); CHKERRQ(ierr); \
    ierr = PetscOptionsHasName(NULL, NULL, "-D2", &check_arg_D); CHKERRQ(ierr); \
    if (check_arg_J && check_arg_D) {					\
      ierr = PetscOptionsGetReal(NULL, NULL, "-J2", &J2, &check_arg_J); CHKERRQ(ierr); \
      ierr = PetscOptionsGetReal(NULL, NULL, "-D2", &Delta2, &check_arg_D); CHKERRQ(ierr); \
    } else {								\
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the value of the next-nearest neighbours interaction with option -J2 and anisotropy with option -D2\n"); CHKERRQ(ierr); \
      return(1);							\
    }									\
  }
#else
#define parse_args()							\
  PetscReal J1 = 0.0;							\
  PetscReal Delta1 = 0.0;						\
  PetscReal J2 = 0.0;							\
  PetscReal Delta2 = 0.0;						\
  PetscBool check_arg_J = PETSC_FALSE,					\
    check_arg_D = PETSC_FALSE,						\
    check_arg_nn = PETSC_FALSE,						\
    check_arg_nnn = PETSC_FALSE;					\
									\
  ierr = PetscOptionsHasName(NULL, NULL, "-nn", &check_arg_nn); CHKERRQ(ierr); \
  ierr = PetscOptionsHasName(NULL, NULL, "-nnn", &check_arg_nnn); CHKERRQ(ierr); \
  									\
  if (!check_arg_nn && !check_arg_nnn) {				\
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide at least one kind of interaction between spins:\n -nn: nearest neighbours\n -nnn: next-nearest neighbours\n"); CHKERRQ(ierr); \
    return(1);								\
  }									\
									\
  /*									\
  PetscInt  _latT, Lx, Ly;						\
  PetscBool check_arg_lat,						\
  check_arg_lat_x,							\
  check_arg_lat_y;							\
  ierr = PetscOptionsHasName(NULL, NULL, "-lat", &check_arg_lat); CHKERRQ(ierr); \
  if (!check_arg_lat) {							\
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the lattice type with option -lat:\n <0> chain1D, <1> square2D \n"); \
    return 1;								\
  }									\
    if (check_arg_lat) {						\
      ierr = PetscOptionsGetInt(NULL, NULL, "-lat", &_latT, &check_arg_lat); CHKERRQ(ierr); \
      if (_latT < 0 || _latT > 1) {					\
	ierr = PetscPrintf(PETSC_COMM_WORLD, "The allowed options are:\n <0> chain1D, <1> square2D \n"); \
	return 1;							\
      }									\
      if (_latT == 1) {							\
	ierr = PetscOptionsHasName(NULL, NULL, "-Lx", &check_arg_lat_x); CHKERRQ(ierr); \
	ierr = PetscOptionsHasName(NULL, NULL, "-Ly", &check_arg_lat_y); CHKERRQ(ierr); \
	if (check_arg_lat_x && check_arg_lat_y) {			\
	  ierr = PetscOptionsGetInt(NULL, NULL, "-Lx", &Lx, &check_arg_lat_x); CHKERRQ(ierr); \
	  ierr = PetscOptionsGetInt(NULL, NULL, "-Ly", &Ly, &check_arg_lat_y); CHKERRQ(ierr); \
	}								\
	else {								\
	  ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the dimensions of the square lattice with -Lx and -Ly\n"); CHKERRQ(ierr); \
	  return 1;							\
	}								\
      }									\
    }									\
*/									\
									\
  if (check_arg_nn) {							\
    ierr = PetscOptionsHasName(NULL, NULL, "-J1", &check_arg_J); CHKERRQ(ierr); \
    ierr = PetscOptionsHasName(NULL, NULL, "-D1", &check_arg_D); CHKERRQ(ierr); \
    if (check_arg_J && check_arg_D) {					\
      ierr = PetscOptionsGetReal(NULL, NULL, "-J1", &J1, &check_arg_J); CHKERRQ(ierr); \
      ierr = PetscOptionsGetReal(NULL, NULL, "-D1", &Delta1, &check_arg_D); CHKERRQ(ierr); \
    } else {								\
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the value of the nearest neighbours interaction with option -J1 and anisotropy with option -D1\n"); CHKERRQ(ierr); \
      return(1);							\
    }									\
  }									\
  									\
  if (check_arg_nnn) {							\
    ierr = PetscOptionsHasName(NULL, NULL, "-J2", &check_arg_J); CHKERRQ(ierr); \
    ierr = PetscOptionsHasName(NULL, NULL, "-D2", &check_arg_D); CHKERRQ(ierr); \
    if (check_arg_J && check_arg_D) {					\
      ierr = PetscOptionsGetReal(NULL, NULL, "-J2", &J2, &check_arg_J); CHKERRQ(ierr); \
      ierr = PetscOptionsGetReal(NULL, NULL, "-D2", &Delta2, &check_arg_D); CHKERRQ(ierr); \
    } else {								\
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Need to provide the value of the next-nearest neighbours interaction with option -J2 and anisotropy with option -D2\n"); CHKERRQ(ierr); \
      return(1);							\
    }									\
  }
#endif

#endif
