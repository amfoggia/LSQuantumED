#include "correlation.hpp"
#include "hamiltonian.hpp"

/* --------------------------------------------------------------------------- */
// ---------------------------- Correlation Class ---------------------------- //
/* --------------------------------------------------------------------------- */

PetscErrorCode Correlation::OpOnStateVector(Environment& env,
					    PetscInt spin_i,
					    PetscInt spin_j,
					    Vec& state,
					    Vec& rhs) {
  
  PetscErrorCode ierr = 0;
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "CorrOpOnStateVector"};
#endif
    
    Mat OpIJ;
    
    ierr = OpCreate(env, spin_i, spin_j, OpIJ);
    ierr = MatMult(OpIJ, state, rhs); CHKERRQ(ierr);
    
    ierr = MatDestroy(&OpIJ); CHKERRQ(ierr);
    
#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("CorrOpOnStateVector");
#endif
  
  return ierr;
}

/* --------------------------------------------------------------------------- */

PetscScalar Correlation::ExpectVal(Environment& env,
				   PetscInt spin_i,
				   PetscInt spin_j,
				   Vec& state) {

  PetscScalar exp_val = 0.0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "CorrExpectVal"};
#endif

    PetscErrorCode ierr = 0;
    Vec rhs, state_conj;

    ierr = VecDuplicate(state, &rhs); CHKERRQ(ierr);
    ierr = VecDuplicate(state, &state_conj); CHKERRQ(ierr);

    ierr = VecCopy(state, state_conj); CHKERRQ(ierr);
    ierr = VecConjugate(state_conj); CHKERRQ(ierr);
    ierr = OpOnStateVector(env, spin_i, spin_j, state, rhs); CHKERRQ(ierr);
    ierr = VecDot(state, rhs, &exp_val); CHKERRQ(ierr);

    ierr = VecDestroy(&state_conj); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs); CHKERRQ(ierr);

#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("CorrExpectVal");
#endif
  
  return exp_val;
}

/* --------------------------------------------------------------------------- */
// -------------------------------- SzSz Class ------------------------------- //
/* --------------------------------------------------------------------------- */

SzSz::SzSz(Environment& m_env,
	   Basis& m_basis)
  : diag_value{0},
    mpi_size{m_env.mpi_size},
    mpi_rank{m_env.mpi_rank},
    basis{m_basis}
{}

/* --------------------------------------------------------------------------- */

void SzSz::OpOnBasisElems(Environment& env,
			  PetscInt spin_i,
			  PetscInt spin_j,
			  PetscInt basis_elem) {

  PetscInt bitA, bitB;
  bitA = (basis_elem >> spin_i) & 0x1;
  bitB = (basis_elem >> spin_j) & 0x1;

  if (bitA ^ bitB)
    diag_value = -1;
  else
    diag_value = 1;
}

/* --------------------------------------------------------------------------- */

PetscErrorCode SzSz::OpOnStateVector(Environment& env,
				     PetscInt spin_i,
				     PetscInt spin_j,
				     Vec& state,
				     Vec& rhs) {
  
  PetscErrorCode ierr = 0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "SzSzOpOnStateVector"};
#endif

    PetscScalar * rhs_array;
    PetscInt local_size, rest;

    rest = basis.size%mpi_size;
    local_size = basis.size/mpi_size + (rest > mpi_rank);
  
    ierr = VecCopy(state, rhs); CHKERRQ(ierr);
    ierr = VecScale(rhs, 0.25); CHKERRQ(ierr);
    ierr = VecGetArray(rhs, &rhs_array); CHKERRQ(ierr);

    for (PetscInt elem = 0; elem < local_size; ++elem) {
      OpOnBasisElems(env, spin_i,spin_j,basis.int_basis[elem]);
      rhs_array[elem] *= diag_value;
    }
  
    ierr = VecRestoreArray(rhs, &rhs_array); CHKERRQ(ierr);

#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("SzSzOpOnStateVector");
#endif
  
  return ierr;
}

/* --------------------------------------------------------------------------- */
// -------------------------------- SpSm Class ------------------------------- //
/* --------------------------------------------------------------------------- */

SpSm::SpSm(Environment& m_env,
	   Basis& m_basis)
  : coup_elem{0},
    swap{PETSC_FALSE},
    basis{m_basis}
{}

/* --------------------------------------------------------------------------- */

void SpSm::OpOnBasisElems(Environment& env,
			  PetscInt spin_i,
			  PetscInt spin_j,
			  PetscInt basis_elem) {

  PetscInt bitA, bitB;
  PetscInt swapped_basis_elem;
  bitA = (basis_elem >> spin_i) & 0x1;
  bitB = (basis_elem >> spin_j) & 0x1;
  swap = PETSC_FALSE;

  if (spin_i == spin_j) {
    if (bitA == 1) {
      coup_elem = HamiltHelper::search_index(basis_elem, basis.nspins);
      swap = PETSC_TRUE;
    }
  }
  else
    if (bitA == 0 && bitB == 1) {
      swapped_basis_elem = basis_elem ^ (1ull << spin_i);
      swapped_basis_elem ^= (1ull << spin_j);
      coup_elem = HamiltHelper::search_index(swapped_basis_elem, basis.nspins);
      swap = PETSC_TRUE;
    }
}

/* --------------------------------------------------------------------------- */

PetscErrorCode SpSm::OpCreate(Environment& env,
			      PetscInt spin_i,
			      PetscInt spin_j,
			      Mat& OpIJ) {

  PetscErrorCode ierr = 0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "SpSmOpCreate"};
#endif
  
    PetscInt Istart, Iend, local_index;

    ierr = MatCreate(PETSC_COMM_WORLD, &OpIJ); CHKERRQ(ierr);
    ierr = MatSetType(OpIJ, MATMPIAIJ); CHKERRQ(ierr);
    ierr = MatSetSizes(OpIJ, PETSC_DECIDE, PETSC_DECIDE, basis.size, basis.size); CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(OpIJ, 1, NULL, 1, NULL); CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(OpIJ, &Istart, &Iend); CHKERRQ(ierr);
    for (PetscInt elem = Istart; elem < Iend; ++elem) {
      local_index = elem - basis.global_start_index;

      OpOnBasisElems(env, spin_i, spin_j, basis.int_basis[local_index]);
      if (swap == PETSC_TRUE) {
	ierr = MatSetValue(OpIJ,
			   elem,
			   coup_elem,
			   1.0,
			   INSERT_VALUES);
	CHKERRQ(ierr);
      }
    }
    ierr = MatAssemblyBegin(OpIJ,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(OpIJ,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("SpSmOpCreate");
#endif
  
  return ierr;
}

/* --------------------------------------------------------------------------- */
// -------------------------------- SmSp Class ------------------------------- //
/* --------------------------------------------------------------------------- */

SmSp::SmSp(Environment& m_env,
	   Basis& m_basis)
  : SpSm{m_env, m_basis}
{}

/* --------------------------------------------------------------------------- */

void SmSp::OpOnBasisElems(Environment& env,
			  PetscInt spin_i,
			  PetscInt spin_j,
			  PetscInt basis_elem) {
  
  PetscInt bitA, bitB;
  PetscInt swapped_basis_elem;
  bitA = (basis_elem >> spin_i) & 0x1;
  bitB = (basis_elem >> spin_j) & 0x1;
  swap = PETSC_FALSE;

  if (spin_i == spin_j) {
    if (bitA == 0) {
      coup_elem = HamiltHelper::search_index(basis_elem, basis.nspins);
      swap = PETSC_TRUE;
    }
  }
  else
    if (bitA == 1 && bitB == 0) {
      swapped_basis_elem = basis_elem ^ (1ull << spin_i);
      swapped_basis_elem ^= (1ull << spin_j);
      coup_elem = HamiltHelper::search_index(swapped_basis_elem, basis.nspins);
      swap = PETSC_TRUE;
    }
}

/* --------------------------------------------------------------------------- */
// -------------------------------- SiSj Class ------------------------------- //
/* --------------------------------------------------------------------------- */

SiSj::SiSj(Environment& m_env,
	   Basis& m_basis)
  : Cpm{SpSm{m_env,m_basis}},
    Cmp{SmSp(m_env,m_basis)},
    Czz{SzSz{m_env,m_basis}}
{}

/* --------------------------------------------------------------------------- */

PetscScalar SiSj::ExpectVal(Environment& env,
			    PetscInt spin_i,
			    PetscInt spin_j,
			    Vec& state) {

  PetscScalar exp_val = 0.0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "SiSjExpectVal"};
#endif
  
    exp_val += Cpm.ExpectVal(env, spin_i,spin_j,state);
    exp_val += Cmp.ExpectVal(env, spin_i,spin_j,state);
    exp_val /= 2.0;
    exp_val += Czz.ExpectVal(env, spin_i,spin_j,state);

#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("SiSjExpectVal");
#endif
  
  return exp_val;
}
