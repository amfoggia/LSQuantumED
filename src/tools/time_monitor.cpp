#include "meson_config.h"
#include "time_monitor.hpp"

/* --------------------------------------------------------------------------- */
/* -------------------------- TimeMonitor Class ------------------------------ */
/* --------------------------------------------------------------------------- */

Tools::TimeMonitor::TimeMonitor(MPI_Comm& m_comm)
  : info{std::map<std::string,tData>{}},
    comm{m_comm} {}

/* --------------------------------------------------------------------------- */

PetscErrorCode Tools::TimeMonitor::PrintTimeInfoFull() {
  PetscErrorCode ierr = 0;
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     "\n===============================================================================\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     "----------------------------- Execution Times ---------------------------------\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\n              Function:    Current ET      Total ET   Calls\n" ); CHKERRQ(ierr);
  for (auto it = info.cbegin(); it != info.cend(); ++it)
    ierr = PetscPrintf(PETSC_COMM_WORLD, "* %20s:   %10e   %10e   %d\n",
		       it->first.c_str(),
		       it->second.current_elapsed_time.count(),
		       it->second.total_elapsed_time.count(),
		       it->second.n_calls); CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     "===============================================================================\n"); CHKERRQ(ierr);
  return(ierr);
}

/* --------------------------------------------------------------------------- */

PetscErrorCode Tools::TimeMonitor::PrintTimeInfoFunc(std::string fn_label) {

  if (info.count(fn_label) == 0)
    throw std::invalid_argument("The function label provided (to print time information) does not exist.");
  
  PetscErrorCode ierr = 0;
  ierr = PetscPrintf(PETSC_COMM_WORLD, "* %20s: %10e -- %10e -- %d\n",
		     fn_label.c_str(),
		     info[fn_label].current_elapsed_time.count(),
		     info[fn_label].total_elapsed_time.count(),
		     info[fn_label].n_calls); CHKERRQ(ierr);
  return(ierr);
}

/* --------------------------------------------------------------------------- */

Tools::ScopedTimer::ScopedTimer(TimeMonitor& m_tm, std::string m_func)
  : start{Clock::now()},
    tm{m_tm},
    func{m_func}
{ tm.info[func].n_calls += 1; }

/* --------------------------------------------------------------------------- */

Tools::ScopedTimer::~ScopedTimer() {
  auto stop = Clock::now();
  tm.info[func].current_elapsed_time = std::chrono::duration<PetscReal>{stop - start};
  tm.info[func].total_elapsed_time += tm.info[func].current_elapsed_time;
}
