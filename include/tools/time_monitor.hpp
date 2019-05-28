#ifndef __TIME_MONITOR_H
#define __TIME_MONITOR_H

#include "meson_config.h"
#include <petscsys.h>
#include <chrono>
#include <map>
#include <iostream>
#include <mpi.h>

/* --------------------------------------------------------------------------- */
// --------------------------- Timing function ------------------------------- //
/* --------------------------------------------------------------------------- */

// Function/Machinery to time the code in a "clean" way
// The idea is that there will be only one monitoring object
// and its data is going to be filled each time the ScopedTimer gets
// called. It time it gets called it initializes
// the std::chrono::steady_clock::now() and when it gets destroyed
// (after going out of scope) it calls the same function again and
// return the elapsed time.


#ifdef TIME_CODE
#define TIME_IT(tm,func_name,func)		\
  {						\
    Tools::ScopedTimer _timer_{tm,func_name};	\
    (func);					\
  }
#else
#define TIME_IT(tm,func_name,func)		\
  (func);
#endif

/* --------------------------------------------------------------------------- */
// --------------------------- Tools namespace ------------------------------- //
/* --------------------------------------------------------------------------- */

/**
 * @namespace Tools
 * @brief Extra functions relevant to the simulation.
 * This includes timing functions and random-generator functions.
 */

namespace Tools {

  /**
   * @struct tData
   * @brief Time data collector.
   */
  struct tData {
    std::chrono::duration<PetscReal> current_elapsed_time; /**< Execution time of one call to a function. */
    std::chrono::duration<PetscReal> total_elapsed_time; /**< Execution time of all the calls to a function. */
    PetscInt n_calls; /**< Number of times a function gets called. */
  };

  /**
   * @class TimeMonitor
   * @brief Useful functions to print time-related data.
   */
  class TimeMonitor {
  
    std::map<std::string,tData> info; /**< Time-related information about a function. */
    MPI_Comm comm; /**< MPI communicator. */

  public:

    /**
     * @fn TimeMonitor(MPI_Comm&)
     * @brief Constructor.
     * @param[in] m_comm MPI communicator.
     */
    TimeMonitor(MPI_Comm& m_comm = PETSC_COMM_WORLD);
    friend class ScopedTimer;

    /**
     * @fn PrintTimeInfoFull
     * @brief Prints time-related information of all functions that were time-tested.
     * @return Error value.
     */
    PetscErrorCode PrintTimeInfoFull();

    /**
     * @fn PrintTimeInfoFunc(std::string)
     * @brief Prints time-related information of a function.
     * @return Error value.
     */
    PetscErrorCode PrintTimeInfoFunc(std::string fn_label);
  };

  /**
   * @class ScopedTimer
   * @brief Initializes and stops a clock for timing the code.
   */
  class ScopedTimer {

    using Clock = std::chrono::steady_clock; /**< Shortcut for std::chrono::steady_clock. */
    std::chrono::time_point<Clock,std::chrono::nanoseconds> start; /**< Start of the clock. */
    TimeMonitor& tm; /**< TimeMonitor object. */
    std::string func; /**< Name of the function to monitor. */
  
  public:

    /**
     * @fn ScopedTimer(TimeMonitor&, std::string)
     * @brief Constructor.
     * @param[in] m_tm TimeMonitor object.
     * @param[in] m_func Name of the function to monitor.
     */
    ScopedTimer(TimeMonitor& m_tm, std::string m_func);
    ~ScopedTimer();
  };

} // -- Namespace

#endif
