#ifndef __RANDOM_DISORDER_H
#define __RANDOM_DISORDER_H

#include "meson_config.h"
#include <petscsys.h>
#include <map>
#include <mpi.h>
#include <random>

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
   * @class RandomDisorder
   * @brief Creates the random disorder constant.
   */
  class RandomDisorder {
    
  public:
    PetscReal min; /**< Minimum value for disorder. */
    PetscReal max; /**< Maximum value for disorder. */
    PetscInt rep; /**< Number of repetitions for disorder computations. */
    PetscInt nspins; /**< Number of spins in the system. */
    PetscInt* srand_seq; /**< Vector with seeds for the random generator. */
    PetscReal* hi; /**< Vector with disorer values per site. */

    /**
     * @fn RandomDisorder(PetscInt, PetscBool, PetscReal, PetscReal, PetscInt)
     * @brief Constructor. Constructs the (random) seeds and the hi disoder vector.
     * @param[in] m_nspins Number of spins in the system.
     * @param[in] m_reprod Should the seeds be random numbers as well or reproducible?
     * @param[in] m_min Minimum value for disorder.
     * @param[in] m_max Maximum value for disorder.
     * @param[in] m_rep Number of repetitions for disorder average.
     */
    RandomDisorder(PetscInt m_nspins,
    		   PetscBool m_reprod,
    		   PetscReal m_min,
    		   PetscReal m_max,
    		   PetscInt m_rep);
    
    /**
     * @fn ~RandomDisorder
     * @brief Destructor. Destroys the PETSc Mat object.
     */
    ~RandomDisorder();
    
    /**
     * @fn srand_seq_init(PetscInt, PetscInt, PetscInt, PetscInt)
     * @brief Creates a list of seeds for the disorer average.
     * @param[in] reprod Should the seeds be random numbers as well or reproducible?
     * @return Vector with seeds.
     */
    PetscInt* srand_seq_init(PetscBool reprod);
    
    /**
     * @fn hi_init(PetscInt)
     * @brief Initializes the hi vector.
     * @param[in] dis_index Index of the srand_seq vector to use as seed.
     * Each element corresponds to the random disorder of each spin.
     */
    void hi_init(PetscInt dis_index);
  };

} // -- Namespace

#endif
