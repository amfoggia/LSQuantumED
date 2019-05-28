#include "environment.hpp"
#include <stdexcept>

/* --------------------------------------------------------------------------- */
// --------------------------- Constructor ----------------------------------- //
/* --------------------------------------------------------------------------- */

Environment::Environment(int argc,
			 char ** argv,
			 PetscInt m_nspins,
			 char* m_help)
{
  if (m_nspins%2 != 0)
    throw std::invalid_argument("The total number of spins has to be EVEN.");

  if (m_nspins <= 0)
    throw std::invalid_argument("The total number of spins has to be POSITIVE.");

  ierr = SlepcInitialize(&argc, &argv, (char*) 0, m_help);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
#ifdef TIME_CODE
  tm = Tools::TimeMonitor{PETSC_COMM_WORLD};
#endif

  nspins = m_nspins;
    
  if (ierr != 0)
    throw std::runtime_error("Slepc was not correctly initialized\n");
}

