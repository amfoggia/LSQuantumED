#ifndef __BASIS_H
#define __BASIS_H

#include "environment.hpp"
#include "meson_config.h"
#include <boost/dynamic_bitset.hpp>
#include <petscsys.h>
#include <iostream>
#include <cmath>
#include <vector>

/**
 * @namespace BasisHelper
 * @brief Helper functions for the Basis class.
 */
namespace BasisHelper {

  /**
   * @fn BasisHelper::search_elem(PetscInt, PetscInt, PetscInt)
   * @brief Given an index of the basis it retireves the corresponding element.
   * @param[in] index Index of the basis.
   * @param[in] nspins Number of spins in the system. 
   * @param[in] nspins_up Number of spins up. 
   * @return Element.
   */
  PetscInt search_elem(PetscInt index,
		       PetscInt nspins,
		       PetscInt nspins_up);
}

/**
 * @class Basis
 * @brief Creates the basis elements for a system with nspins and constant total magnetization.
 */

class Basis {

private:
  PetscInt max_mag; /**< Maximum magnetization allowed for the system. */

  /**
   * @fn compute_size
   * @brief Computes the number of elements in the basis.
   * @return Number of elements of basis.
   */
  PetscInt compute_size(); 

  /**
   * @fn generate_int_basis
   * @brief Creates the 64-bit integer basis elements and stores them in an std::vector. 
   * @return Vector with basis elements.
   */
  std::vector<PetscInt>& generate_int_basis();

#ifdef DEVEL
  /**
   * @fn generate_bit_basis
   * @brief Creates the bit representationof the basis elements and stores them in an std::vector. 
   * @return Vector with basis elements.
   */
  std::vector<boost::dynamic_bitset<>>& generate_bit_basis(); // TO DO: DELETE 
#endif
  
  // ALL THE ELEMENTS SHOULD BE MUTABLE CONST AND AL THE FUNCTIONS (NOT CONSTRUCTORS) SHOULD BE 
  // CONST (AT THE END)
public:
  PetscInt nspins; /**< Number of spins in the system. */
  PetscInt total_mag; /**< Total magnetization in the system. */
  PetscInt nspins_up; /**< Number of spins pointing up in the system. */
  PetscInt size; /**< Number of basis elements. */
  PetscMPIInt mpi_size; /**< Rank label for each process. */
  PetscMPIInt mpi_rank; /**< Number of MPI processes. */
  PetscInt local_size; /**< Number of basis elements in each process. */
  PetscInt global_start_index; /**< Global index of the first basis element of each process. */
  std::vector<PetscInt> int_basis; /**< Basis elements. */
#ifdef DEVEL
  std::vector<boost::dynamic_bitset<>> bit_basis; /**< Bit representation of the basis elements. */ // TO DO: DELETE
#endif
  
  // Constructors
  /**
   * @fn Basis
   * @brief Default constructor.
   */
  Basis();

  /**
   * @fn Basis(Environment&, PetscInt)
   * @brief Constructor.
   * @param[in] env Environment object.
   * @param[in] total_mag Total magnetization of the system.
   */
  Basis(Environment& env,PetscInt total_mag);

  /**
   * @fn Basis(const Basis&)
   * @brief Copy constructor.
   * @param[in] rhs Basis object.
   */
  Basis(const Basis& rhs);

  /**
   * @fn Basis(Basis&&)
   * @brief Move constructor.
   * @param[in] rhs Basis object.
   */
  Basis(Basis&& rhs) noexcept;

  // Assignment operators
  /**
   * @fn operator=(const Basis& rhs)
   * @brief Copy assignment.
   * @param[in] rhs Basis object.
   * @return Basis object.
   */
  Basis& operator=(const Basis& rhs);

  /**
   * @fn operator=(Basis&& rhs)
   * @brief Move assignment.
   * @param[in] rhs Basis object.
   * @return Basis object.
   */
  Basis& operator=(Basis&& rhs) noexcept;
  
};

#endif
