#ifndef __ENVIRONMENT_H
#define __ENVIRONMENT_H

#include "meson_config.h"
#include "time_monitor.hpp"
#include <iostream>
#include <slepcsys.h>

/**
 * @class Environment
 * @brief Initializes and finalizes the SLEPc environment (SlepcInitialize-SlepcFinalize). It sets te number of spins for it to be consistent throughout the code.
 */

class Environment {

private:
  PetscErrorCode ierr = 0; /**< Variable for the return error value of PETSc. */
  
public:
  PetscInt nspins; /**< Number of spins in the system. */
  PetscMPIInt mpi_rank; /**< Rank label for each process. */
  PetscMPIInt mpi_size; /**< Number of MPI processes. */
#ifdef TIME_CODE
  Tools::TimeMonitor tm; /**< TimeMonitor object to measure time across the code. */
#endif
  PetscBool disorder_flg; /**< Flag that states if disorder realizations are going to be performed. */
  
  // Constructors
  /**
   * @fn Environment(int, char**, char[])
   * @brief Constructor.
   * Initializes the SLEPc environment.
   * @param[in] argc Number of arguments given to the executable.
   * @param[in] argv Array of arguments given to the executable.
   * @param[in] help String that describes the purpose of the executable.
   */
  Environment(int argc,
	      char ** argv,
	      char help[] = nullptr);

  /**
   * @fn Environment(int, char**, PetscInt, char[])
   * @brief Constructor.
   * Initializes the SLEPc environment.
   * @param[in] argc Number of arguments given to the executable.
   * @param[in] argv Array of arguments given to the executable.
   * @param[in] m_nspins Number of spins in the system.
   * @param[in] help String that describes the purpose of the executable.
   */
  Environment(int argc,
	      char ** argv,
	      PetscInt m_nspins,
	      char help[] = nullptr);

  // Destructor
  /**
   * @fn ~Environment
   * @brief Destructor.
   * Finalizes the SLEPc environment.
   */
  ~Environment() { ierr = SlepcFinalize(); }
  
};

#endif
