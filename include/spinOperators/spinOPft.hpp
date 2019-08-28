#ifndef __SPINOPFT_H
#define __SPINOPFT_H

#include "basis.hpp"
#include "lattice.hpp"
#include "hamiltonian.hpp"

/**
 * @struct Sq_data
 * @brief Necessary data for the construction of the Sq operator.
 * Later, in addition of other things, needed for the computation of the Dynamical Structure Factor.
 */
struct Sq_data {
  Basis* b0; /**< Basis of the Mag = 0 subspace. */
  Basis* b1; /**< Basis of the new subspace. */
};

/* --------------------------------------------------------------------------- */

/**
 * @class Sq
 * @brief Creates the Fourier transform of the spin operator.
 * Provides functions to apply the operator to a state vector and to a basis element.
 */

class Sq {
  
public:

  /**
   * @fn OpOnBasisElems(Environment&, PetscInt, PetscInt, PetscInt, PetscInt)
   * @brief Computes the operators applied on the basis elements.
   * Each particular correlation (derived classes) will compute it in a different way.
   * This is what a user should provide if it were to add a new correlation function.
   * @param[in] env Environment object.
   * @param[inout] basis_elem Element of basis to which one applies the operator.
   * @param[in] spin Lattice site to which the operator is applied.
   */
  virtual void OpOnBasisElems(Environment& env, PetscInt basis_elem, PetscInt spin) = 0;

  /**
   * @fn OpOnStateVector(Environment&, Vec&, Vec&, PetscInt)
   * @brief Computes the operator applied to a state vector (linear combination of basis elements).
   * @param[in] env Environment object.
   * @param[in] state State vector to which one applies the operator.
   * @param[inout] rhs Vector obtained after applying the operator on the right hand side vector.
   * @param[in] spin Lattice site to which the operator is applied.
   * @return Error value.
   */
  virtual PetscErrorCode OpOnStateVector(Environment& env, Vec& state, Vec& rhs, PetscInt spin) = 0;

   /**
   * @fn ~Sq
   * @brief Virtual destructor.
   */
  virtual ~Sq() {}
};

/* --------------------------------------------------------------------------- */

/**
 * @class Sqz
 * @brief Creates the Fourier transform of the spin operator Sz.
 * Provides functions to apply the operator to a state vector and to a basis element.
 */

class Sqz : public Sq {

private:
  PetscInt nspins; /**< Number of spins in the system. */
  Sq_data& data; /**< Sq_data object. */
  PetscMPIInt mpi_size; /**< Number of MPI processes. */
  PetscMPIInt mpi_rank; /**< Rank label for each process. */

#ifdef DEVEL
public:
#endif
  PetscReal value; /**< Coefficient obtained after applying the operator to a basis element. */
  
public:

  /**
   * @fn Sqz(Environment&)
   * @brief Constructor.
   * @param[in] m_env Environment object.
   * @param[in] m_data Sq_data object.
   */
  Sqz(Environment& m_env, Sq_data& m_data);

  void OpOnBasisElems(Environment& env, PetscInt basis_elem, PetscInt spin);

  PetscErrorCode OpOnStateVector(Environment& env, Vec& state, Vec& rhs, PetscInt spin);
 
  ~Sqz() {}
};

/* --------------------------------------------------------------------------- */

/**
 * @class Sqp
 * @brief Creates the Fourier transform of the spin operator Sp.
 * Provides functions to apply the operator to a state vector and to a basis element.
 */

class Sqp : public Sq {

private:
  PetscInt nspins; /**< Number of spins in the system. */
  Sq_data& data; /**< Sq_data object. */
  PetscMPIInt mpi_size; /**< Number of MPI processes. */
  PetscMPIInt mpi_rank; /**< Rank label for each process. */

#ifdef DEVEL
public:
#endif
  PetscInt value; /**< Basis element to which a particular elements is coupled after applying the operator Sp. */
  
public:

  /**
   * @fn Sqp(Environment&)
   * @brief Constructor.
   * @param[in] m_env Environment object.
   * @param[in] m_data Sq_data object.
   */
  Sqp(Environment& m_env, Sq_data& m_data);

  void OpOnBasisElems(Environment& env, PetscInt basis_elem, PetscInt spin);

  PetscErrorCode OpOnStateVector(Environment& env, Vec& state, Vec& rhs, PetscInt spin);

  ~Sqp() {}
};

#include "spinOPft.tmpl_func.hpp"

#endif
