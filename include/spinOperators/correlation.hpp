#ifndef __CORRELATIONS_H
#define __CORRELATIONS_H

#include "basis.hpp"

/**
 * @class Correlation
 * @brief Pure virtual class for the creation of the correlation functions.
 */

class Correlation {

public:
  
  /**
   * @fn OpOnBasisElems(Environment&, PetscInt, PetscInt, PetscInt)
   * @brief Computes the operators applied on the basis elements.
   * Each particular correlation (derived classes) will compute it in a different way.
   * This is what a user should provide if it were to add a new correlation function.
   * @param[in] env Environment object.
   * @param[in] spin_i First spin involved in the correlation function.
   * @param[in] spin_j Second spin involved in the correlation function.
   * @param[inout] basis_elem Element of basis to which one applies the operator.
   */
  virtual void OpOnBasisElems(Environment& env, PetscInt spin_i, PetscInt spin_j, PetscInt basis_elem) = 0;

  /**
   * @fn OpCreate(Environment&, PetscInt, PetscInt)
   * @brief Creates the matrix (Mat object) corresponding to the operator.
   * @param[in] env Environment object.
   * @param[in] spin_i First spin involved in the correlation function.
   * @param[in] spin_j Second spin involved in the correlation function.
   * @param[inout] OpIJ Matrix representing the operator.
   * @return Error value.
   */
  virtual PetscErrorCode OpCreate(Environment& env, PetscInt spin_i, PetscInt spin_j, Mat& OpIJ) = 0;

  /**
   * @fn OpOnStateVector(Environment&, PetscInt, PetscInt, Vec&, Vec&)
   * @brief Computes the operator applied to a state vector (linear combination of basis elements).
   * A default implementation is given, which includes creating the matrix of the operator and
   * using Petsc functions for do the matrix-vector product.
   * @param[in] env Environment object.
   * @param[in] spin_i First spin involved in the correlation function.
   * @param[in] spin_j Second spin involved in the correlation function.
   * @param[in] state State vector to use in the computation of the expectation value of the operator.
   * @param[inout] rhs Vector obtained after applying the operator on the right hand side vector.
   * @return Error value.
   */
  virtual PetscErrorCode OpOnStateVector(Environment& env, PetscInt spin_i, PetscInt spin_j, Vec& state, Vec& rhs);

  /**
   * @fn ExpectVal(Environment&, PetscInt, PetscInt, Vec&)
   * @brief Computes the expectation value of the operator.
   * A default implementation is given. It uses a Petsc function to 
   * compute the vector-vector dot product.
   * @param[in] env Environment object.
   * @param[in] spin_i First spin involved in the correlation function.
   * @param[in] spin_j Second spin involved in the correlation function.
   * @param[in] state State vector to use in the computation of the expectation value of the operator.
   * @return Expectation value. 
   */
  virtual PetscScalar ExpectVal(Environment& env, PetscInt spin_i, PetscInt spin_j, Vec& state);

  /**
   * @fn ~Correlation
   * @brief Virtual destructor.
   */
  virtual ~Correlation() {}
};

/* --------------------------------------------------------------------------- */

/**
 * @class SzSz
 * @brief Computes the SzSz correlation function of two spins.
 */

class SzSz : public Correlation {

private:
  PetscInt diag_value; /**< Sign that results from applying the SzSz operator to a basis element. */
  PetscMPIInt mpi_size; /**< Number of MPI processes. */
  PetscMPIInt mpi_rank; /**< Rank label for each process. */
  Basis basis; /**< Basis object. */
  
public:
  
  /**
   * @fn SzSz(Environment&, Basis&)
   * @brief Constructor.
   * @param[in] m_env Environment object.
   * @param[in] m_basis Basis object.
   */
  SzSz(Environment& m_env, Basis& m_basis);
  
  void OpOnBasisElems(Environment& env, PetscInt spin_i, PetscInt spin_j, PetscInt basis_elem);

  PetscErrorCode OpCreate(Environment& env, PetscInt spin_i, PetscInt spin_j, Mat& OpIJ) { return 0; }

  PetscErrorCode OpOnStateVector(Environment& env, PetscInt spin_i, PetscInt spin_j, Vec& state, Vec& rhs);

  ~SzSz() {}
};

/* --------------------------------------------------------------------------- */

/**
 * @class SpSm
 * @brief Computes the S+S- correlation function of two spins.
 */

class SpSm : public Correlation {

protected:
  PetscInt coup_elem; /**< Basis element that results from applying the SpSm operator to another basis element. */
  PetscBool swap; /**< Has it been a swap or not?. */
  Basis basis; /**< Basis object. */
  
public:
  
  /**
   * @fn SpSm(Environment&, Basis&)
   * @brief Constructor.
   * @param[in] m_env Environment object.
   * @param[in] m_basis Basis object.
   */
  SpSm(Environment& m_env, Basis& m_basis);
  
  virtual void OpOnBasisElems(Environment& env, PetscInt spin_i, PetscInt spin_j, PetscInt basis_elem);

  virtual PetscErrorCode OpCreate(Environment& env, PetscInt spin_i, PetscInt spin_j, Mat& OpIJ);

  virtual ~SpSm() {}
};

/* --------------------------------------------------------------------------- */

/**
 * @class SmSp
 * @brief Computes the S-S+ correlation function of two spins.
 */

class SmSp : public SpSm {

public:

  /**
   * @fn SmSp(Environment&, Basis&)
   * @brief Constructor.
   * @param[in] m_env Environment object.
   * @param[in] m_basis Basis object.
   */
  SmSp(Environment& m_env, Basis& m_basis);
  
  void OpOnBasisElems(Environment& env, PetscInt spin_i, PetscInt spin_j, PetscInt basis_elem);

  ~SmSp() {}
};

/* --------------------------------------------------------------------------- */

/**
 * @class SiSj
 * @brief Computes the (Si * Sj) correlation function of two spins.
 */

class SiSj : public Correlation{

private:
  SpSm Cpm; /**< SpSm object. */
  SmSp Cmp; /**< SmSp object. */
  SzSz Czz; /**< SzSz object. */
  
public:

  /**
   * @fn SiSj(Environment&, Basis&)
   * @brief Constructor.
   * @param[in] m_env Environment object.
   * @param[in] m_basis Basis object.
   */
  SiSj(Environment& m_env, Basis& m_basis);

  void OpOnBasisElems(Environment& env, PetscInt spin_i, PetscInt spin_j, PetscInt basis_elem) {}

  PetscErrorCode OpCreate(Environment& env, PetscInt spin_i, PetscInt spin_j, Mat& OpIJ) {return 0;}

  PetscScalar ExpectVal(Environment& env, PetscInt spin_i, PetscInt spin_j, Vec& state);

  ~SiSj() {}
};

#endif
