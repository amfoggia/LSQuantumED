#ifndef __HAMILTONIAN_H
#define __HAMILTONIAN_H

#include "environment.hpp"
#include "meson_config.h"
#include "lattice.hpp"
#include "basis.hpp"
#include "boost/math/special_functions/binomial.hpp"
#include <slepceps.h>

/**
 * @namespace HamiltHelper
 * @brief Helper functions for the Hamiltonian class.
 */

namespace HamiltHelper {
  
  /**
   * @fn HamiltHelper::get_elem_diag_AUX(PetscInt, PetscInt, PetscInt, PetscInt&)
   * @brief Takes two bits, compares them and adds/subtracs according to the result.
   * @param[in] basis_elem Element of the Basis being treated currently in the loop.
   * @param[in] i Current bit.
   * @param[in] spin Bit correspondig to i's neighbour.
   * @param[inout] diag_elem Element of the diagonal.
   */
  void get_elem_diag_AUX(PetscInt basis_elem,
			 PetscInt i,
			 PetscInt spin,
			 PetscInt& diag_elem);

  /**
   * @fn HamiltHelper::get_coup_elems_AUX(Basis&, PetscInt, PetscInt, PetscInt, PetscInt&, PetscInt*)
   * @brief Takes two bits, swaps them and retrieves the basis element corresponding to the new bit configuration.
   * @param[in] basis A Basis object.
   * @param[in] elem Index of the basis element.
   * @param[in] i Current bit.
   * @param[in] spin Bit correspondig to i's neighbour.
   * @param[inout] ncol Number of elements (up to the moment) present in one matrix row (minus the diagonal element).
   * @param[inout] coup_elems Indices of basis elements related to the elem element. 
   */
  void get_coup_elems_AUX(Basis& basis,
			  PetscInt elem,
			  PetscInt i,
			  PetscInt spin,
			  PetscInt& ncol,
			  PetscInt* coup_elems);

  /**
   * @fn HamiltHelper::get_disorder_elem(PetscInt, PetscInt, PetscReal*)
   * @brief Gets diagonal elements of the Hamiltonian matrix corresponding to the disorder part.
   * @param[in] basis_elem Element of the Basis being treated currently in the loop.
   * @param[in] nspins Number of spins.
   * @param[in] hi Random numbers for disorder. The lengh of the vector is nspins.
   * @return Element of the diagonal corresponding to the disorder part.
   */
  PetscReal get_disorder_elem(PetscInt basis_elem,
			      PetscInt nspins,
			      PetscReal* hi);
  
  /**
   * @fn HamiltHelper::get_elem_diag(neigh_type, PetscInt, PetscInt, L&)
   * @brief Gets diagonal elements of the Hamiltonian matrix.
   * @tparam L Lattice type.
   * @param[in] nt Type of neighbour interaction.
   * @param[in] basis_elem Element of the Basis being treated currently in the loop.
   * @param[in] nspins Number of spins.
   * @param[in] lattice Lattice object.
   * @return Element of the diagonal.
   */
  template<typename L>
  PetscInt get_elem_diag(neigh_type nt,
			 PetscInt basis_elem,
			 PetscInt nspins, 
			 L& lattice);
  
  /**
   * @fn HamiltHelper::get_coup_elems(neigh_type, Basis&, PetscInt, L&, PetscInt&, PetscInt*)
   * @brief Gets off-diagonal elements of the Hamiltonian matrix.
   * @tparam L Lattice type.
   * @param[in] nt Type of neighbour interaction.
   * @param[in] basis Basis object.
   * @param[in] elem Index of a basis element.
   * @param[in] lattice Lattice object.
   * @param[in,out] ncol Number of elements present in one matrix row (minus the diagonal element). 
   * @param[in,out] coup_elems Indices of basis elements related to the elem element. 
   */
  template<typename L>
  void get_coup_elems(neigh_type nt,
		      Basis& basis,
		      PetscInt elem,
		      L& lattice,
		      PetscInt& ncol,
		      PetscInt* coup_elems);
  
  /**
   * @fn HamiltHelper::get_elems(neigh_type, Basis&, PetscInt, L&, PetscInt&, PetscInt&, PetscInt*)
   * @brief Gets elements of the Hamiltonian matrix.
   * @tparam L Lattice type.
   * @param[in] nt Type of neighbour interaction.
   * @param[in] basis Basis object.
   * @param[in] elem Index of a basis element.
   * @param[in] lattice Lattice object.
   * @param[in,out] ncol Number of elements present in one matrix row (minus the diagonal element). 
   * @param[in,out] diag_elem Integer value of the diagonal element for one row of the matrix.
   * @param[in,out] coup_elems Indices of basis elements related to the elem element. 
   */
  template<typename L>
  void get_elems(neigh_type nt,
		 Basis& basis,
		 PetscInt elem,
		 L& lattice,
		 PetscInt& ncol,
		 PetscInt& diag_elem,
		 PetscInt* coup_elems);
  
  /**
   * @fn HamiltHelper::search_index(PetscInt, PetscInt)
   * @brief Given a basis element it retrieves its vector index.
   * @param[in] basis_elem Basis element.
   * @param[in] nspins Number of spins in the system. 
   * @return Index.
   */
  PetscInt search_index(PetscInt basis_elem,
			PetscInt nspins);

  /**
   * @fn HamiltHelper::prealloc_info(PetscReal, PetscInt, PetscInt, PetscInt, PetscInt*, PetscInt, PetscInt, PetscInt*, PetscInt*)
   * @brief Gets the vectors with the matrix shape (needed for a PETSc function).
   * @param[in] Delta Value of the anisotropy.
   * @param[in] elem Index of a basis element.
   * @param[in] global_start_index Global index of the first basis element of each process.
   * @param[in] global_end_index Global index of the following to the last basis element of each process.
   * @param[in] coup_elems Indices of basis elements related to the elem element.
   * @param[in] ncol Number of elements present in one matrix row (minus the diagonal element). 
   * @param[in] elem_diag Integer value of the diagonal element for one row of the matrix.
   * @param[in,out] d_nnz Array containing the number of nonzeros in the various rows of the diagonal portion of the local submatrix.
   * @param[in,out] o_nnz Array containing the number of nonzeros in the various rows of the off-diagonal portion of the local submatrix.
   */
  void prealloc_info(PetscReal Delta,
		     PetscInt elem,
		     PetscInt global_start_index,
		     PetscInt global_end_index,
		     PetscInt* coup_elems,
		     PetscInt ncol,
		     PetscInt elem_diag,
		     PetscInt* d_nnz,
		     PetscInt* o_nnz);

}

/**
 * @class Hamiltonian.
 * @brief Construct and holds the Hamiltonian matrix of the system.
 * @tparam L Lattice type.
 */

template<typename L>
class Hamiltonian {

private:
  Basis basis; /**< Basis of the system. */
  L& lattice; /**< Lattice of the system. */
  PetscInt size; /**< Number of basis elements. */
  PetscReal J1; /**< Strength of the nearest neighbour bonds. */
  PetscReal Delta1; /**< Strength of the anisotropy for nearest neighbours. */
  PetscReal J2; /**< Strength of the next-nearest neighbour bonds. */
  PetscReal Delta2; /**< Strength of the anisotropy for next-nearest neighbours. */

  /** 
   * @fn MatInit(Environment&)
   * @brief Set of calls to PETSc function to create the Mat object.
   * @param[in] env Environment object.
   * @return Error value.
   */
  PetscErrorCode MatInit(Environment& env);

  /** 
   * @fn MatClean
   * @brief Call to the PETSc function MatDestroy().
   * @return Error value.
   */
  PetscErrorCode MatClean();
  
public:
  PetscMPIInt mpi_size; /**< Rank label for each process. */
  PetscMPIInt mpi_rank; /**< Number of MPI processes. */
  Mat hamilt;  /**< PETSc matrix object that contains the Hamiltonian of the system. */
  
  // Constructors
  /**
   * @fn Hamiltonian(Environment&, Basis&, L&, PetscReal, PetscReal, PetscReal, PetscReal)
   * @brief Constructor. Constructs the PETSc Mat object.
   * @param[in] env Environment object.
   * @param[in] m_basis Basis object.
   * @param[in] m_lattice Lattice object.
   * @param[in] m_J1 Strenght of the nearest neighbour bonds.
   * @param[in] m_Delta1 Strenght of the nearest neighbours anisotropy.
   * @param[in] m_J2 Strenght of the next-nearest neighbour bonds.
   * @param[in] m_Delta2 Strenght of the next-nearest neighbours anisotropy.
   */
  Hamiltonian(Environment& env,
	      Basis& m_basis,
	      L& m_lattice,
	      PetscReal m_J1,
	      PetscReal m_Delta1,
	      PetscReal m_J2,
	      PetscReal m_Delta2);

  /**
   * @fn ~Hamiltonian
   * @brief Destructor. Destroys the PETSc Mat object.
   */
  ~Hamiltonian() { MatClean(); }

  // Building functions

  /**
   * @fn build_disorder
   * @brief Fills in the diagonal part of the matrix that corresponds to the disorder term in the Hamiltonian.
   * @param[in] env Environment object.
   * @param[in] hi Vector with random numbers simulation disorder.
   * @return Error value.
   */
  PetscErrorCode build_disorder(Environment& env,
				PetscReal* hi);
  
  /**
   * @fn build_diag(Environment&)
   * @brief Fills in the diagonal part of the matrix.
   * @param[in] env Environment object.
   * @return Error value.
   */
  PetscErrorCode build_diag(Environment& env);

  /**
   * @fn build_off_diag(Environment&)
   * @brief Fills in the off-diagonal part of the matrix.
   * @param[in] env Environment object.
   * @return Error value.
   */
  PetscErrorCode build_off_diag(Environment& env);
  
  /**
   * @fn build_hamilt(Environment&)
   * @brief Fills in the matrix.
   * @param[in] env Environment object.
   * @return Error value.
   */
  PetscErrorCode build_hamilt(Environment& env);

  // Retrieve-attributes functions
  /**
   * @fn hamilt_size
   * @brief Gets the linear size of the Hamiltonian. It is the same as the basis number of elements.
   * @return Linear size of the Hamiltonian.
   */
  PetscInt hamilt_size() const { return size; }

  /**
   * @fn hamilt_dim
   * @brief Gets the theoretical size of the Hamiltonian. It is a sparse matrix so there are (way) less elements.
   * @return Theoretical size of the Hamiltonian.
   */
  PetscInt hamilt_dim() const { return size*size; }

  /**
   * @fn hamilt_J1
   * @brief Gets the strength of the nearest neighbour bonds.
   * @return Strength of the nearest neighbour bonds.
   */
  PetscReal hamilt_J1() const { return J1; }

  /**
   * @fn hamilt_J2
   * @brief Gets the strength of the next-nearest neighbour bonds.
   * @return Strength of the nearest next-neighbour bonds.
   */
  PetscReal hamilt_J2() const { return J2; }

  /**
   * @fn hamilt_D1
   * @brief Gets the strength of the anisotropy for nearest neighbours.
   * @return Strength of the anisotropy for nearest neighbours.
   */
  PetscReal hamilt_D1() const { return Delta1; }
  
  /**
   * @fn hamilt_D2
   * @brief Gets the strength of the anisotropy for next-nearest neighbours.
   * @return Strength of the anisotropy for next-nearest neighbours.
   */
  PetscReal hamilt_D2() const { return Delta2; }
  
  /**
   * @fn hamilt_type
   * @brief Gets the type of the Lattice.
   * @return Lattice type.
   */
  lattice_type hamilt_type() const { return lattice.get_type(); }

  // Preallocation
  /**
   * @fn prealloc(Environment&)
   * @brief Calls the PETSc function to preallocate space for the matrix.
   * As it is a sparse matrix, one sets the linear dimentions but when the Mat object os created it does not allocate all that memory, because it is not 
   * going to be needed. If after the creation no preallocation is done, the Mat object is treated as a dynamic array, making the filling process astonishingly slow.\
   * @param[in] env Environment object.
   * @return Error value.
   */
  PetscErrorCode prealloc(Environment& env);
};

#include "hamiltonian.tmpl_func.hpp"

#endif
