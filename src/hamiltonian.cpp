#include "hamiltonian.hpp"
#include "boost/math/special_functions/binomial.hpp"
#include <omp.h>

/* --------------------------------------------------------------------------- */
// ---------------------------- Extra Functions ------------------------------ //
/* --------------------------------------------------------------------------- */

PetscInt HamiltHelper::search_index(PetscInt basis_elem,
				    PetscInt nspins) {

  PetscInt spinup_index, tmp;
  PetscReal index;
  
  index = 0;
  spinup_index = 0;
  tmp = 0;
  for (int bit = 0; bit < nspins; ++bit) {
    tmp = basis_elem >> bit;
    if (tmp & 1) {
      ++spinup_index;
      if (!(bit < spinup_index))
	index += boost::math::binomial_coefficient<double>((double)bit,(double)spinup_index); 
    }
  }
  return (PetscInt)std::floor(index + 0.5);
}

/* --------------------------------------------------------------------------- */

PetscReal HamiltHelper::get_disorder_elem(PetscInt basis_elem,
					  PetscInt nspins,
					  PetscReal* hi) {

  PetscInt bitA;
  PetscReal disorder_elem = 0;

  for (PetscInt i = 0; i < nspins; ++i) { 

    bitA = (basis_elem >> i) & 0x1;

    if (bitA) 
      disorder_elem += hi[i];
    else
      disorder_elem -= hi[i];
  }

  return disorder_elem;
}

/* --------------------------------------------------------------------------- */

void HamiltHelper::get_elem_diag_AUX(PetscInt basis_elem,
				     PetscInt i,
				     PetscInt spin,
				     PetscInt& diag_elem) {

  PetscInt bitA, bitB;

  // Get bits to compare
  bitA = (basis_elem >> i) & 0x1;
  bitB = (basis_elem >> spin) & 0x1;
  
  if (bitA ^ bitB)
    diag_elem -= 1; // If bits are different
  else
    diag_elem += 1; // If bits are equal   
}

/* --------------------------------------------------------------------------- */

void HamiltHelper::get_coup_elems_AUX(Basis& basis,
				      PetscInt elem,
				      PetscInt i,
				      PetscInt spin,
				      PetscInt& ncol,
				      PetscInt* coup_elems) {

  PetscInt bitA, bitB;
  PetscInt basis_elem, swapped_basis_elem;
  basis_elem = basis.int_basis[elem];
  std::vector<PetscInt> swapped_elems;

  // Get the bits to compare
  bitA = (basis_elem >> i) & 0x1;
  bitB = (basis_elem >> spin) & 0x1;
  
  // Compare bits
  if (bitA ^ bitB) {
    
    // Swap bits
    swapped_basis_elem = basis_elem ^ (1ull << i);
    swapped_basis_elem ^= (1ull << spin);
    
    // Search for the element
    coup_elems[ncol] = search_index(swapped_basis_elem, basis.nspins);
    ncol += 1;
  }
}

/* --------------------------------------------------------------------------- */
// ---------------------------- Prealloc Functions --------------------------- //
/* --------------------------------------------------------------------------- */

void HamiltHelper::prealloc_info(PetscReal Delta,
				 PetscInt elem,
				 PetscInt global_start_index,
				 PetscInt global_end_index,
				 PetscInt* coup_elems,
				 PetscInt ncol,
				 PetscInt elem_diag,
				 PetscInt* d_nnz,
				 PetscInt* o_nnz) {

  PetscInt local_index = elem - global_start_index;

  // Diagonal elements
  // if (elem_diag != 0 && std::fabs(Delta) >= 1e-14)
  // d_nnz[local_index] += 1;

  // Off diagonal elements
  for (int i = 0; i < ncol; ++i) {
    if (coup_elems[i] >= global_start_index && coup_elems[i] < global_end_index)
      d_nnz[local_index] += 1;
    else
      o_nnz[local_index] += 1;
  }  
}
