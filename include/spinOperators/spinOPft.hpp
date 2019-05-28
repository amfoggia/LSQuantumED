#ifndef __SPINOPFT_H
#define __SPINOPFT_H

#include "basis.hpp"
#include "lattice.hpp"

/**
 * @struct Sq_data
 * @brief Necessary data for the construction of the Sq operator.
 * Later, in addition of other things, needed for the computation of the Dynamical Structure Factor.
 * @tparam Lattice type.
 */
template<typename L>
struct Sq_data {
  Basis* b0; /**< Basis of the Mag = 0 subspace. */
  Basis* b1; /**< Basis of the new subspace. */
  L* lat; /**< Lattice object. */
  // TODO: qi should be templated on the dimension of the lattice (1D, 2D, 3D)
  std::array<PetscReal,2> qi; /**< Value of q to use. */
};

/* --------------------------------------------------------------------------- */

/**
 * @class Sq
 * @brief Creates the Fourier transform of the spin operator.
 * Provides functions to apply the operator to a state vector and to a basis element.
 * @tparam Lattice type.
 */

template<typename L>
class Sq {
  
public:

  /**
   * @fn OpOnBasisElems(Environment&, PetscInt, PetscInt, PetscInt)
   * @brief Computes the operators applied on the basis elements.
   * Each particular correlation (derived classes) will compute it in a different way.
   * This is what a user should provide if it were to add a new correlation function.
   * @param[in] env Environment object.
   * @param[inout] basis_elem Element of basis to which one applies the operator.
   */
  virtual void OpOnBasisElems(Environment& env, PetscInt basis_elem) = 0;

  /**
   * @fn OpOnStateVector(Environment&, Vec&, Vec&)
   * @brief Computes the operator applied to a state vector (linear combination of basis elements).
   * @param[in] env Environment object.
   * @param[in] state State vector to which one applies the operator.
   * @param[inout] rhs Vector obtained after applying the operator on the right hand side vector.
   * @return Error value.
   */
  virtual PetscErrorCode OpOnStateVector(Environment& env, Vec& state, Vec& rhs) = 0;

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
 * @tparam Lattice type.
 */

template<typename L>
class Sqz : public Sq<L> {

private:
  PetscInt nspins; /**< Number of spins in the system. */
  Sq_data<L>& data; /**< Sq_data object. */
  PetscMPIInt mpi_size; /**< Number of MPI processes. */
  PetscMPIInt mpi_rank; /**< Rank label for each process. */
  std::vector<PetscComplex> exp; /**< Vector with the exponential part of the operator. */

  /**
   * @fn compute_exp
   * @brief Computes the exponential part of the Sqz operator.
   * @return Vector with the exponential part.
   */
  std::vector<PetscComplex>& compute_exp();

#ifdef DEVEL
public:
#endif
  PetscScalar value; /**< Coefficient obtained after applying the operator to a basis element. */
  
public:

  /**
   * @fn Sqz(Environment&)
   * @brief Constructor.
   * @param[in] m_env Environment object.
   * @param[in] m_data Sq_data object.
   */
  Sqz(Environment& m_env, Sq_data<L>& m_data);

  void OpOnBasisElems(Environment& env, PetscInt basis_elem);

  PetscErrorCode OpOnStateVector(Environment& env, Vec& state, Vec& rhs);

  ~Sqz() {}
};

#include "spinOPft.tmpl_func.hpp"

#endif
