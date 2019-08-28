/* --------------------------------------------------------------------------- */
// --------------------------------- Sqz Class ------------------------------- //
/* --------------------------------------------------------------------------- */

Sqz::Sqz(Environment& m_env,
	 Sq_data& m_data)
  : nspins{m_env.nspins},
    data{m_data},
    mpi_size{m_env.mpi_size},
    mpi_rank{m_env.mpi_rank},
    value{0.0}
{}

/* --------------------------------------------------------------------------- */

void Sqz::OpOnBasisElems(Environment& env, PetscInt basis_elem, PetscInt spin) {
  
  PetscInt bitA;
  value = 1.0;
  
  bitA = (basis_elem >> spin) & 0x1;
  if (!bitA)
    value *= -1.0;
}

/* --------------------------------------------------------------------------- */

PetscErrorCode Sqz::OpOnStateVector(Environment& env, Vec& state, Vec& rhs, PetscInt spin) {

  PetscErrorCode ierr = 0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "SqzOpOnStateVector"};
#endif
    
    PetscScalar * rhs_array;
    PetscInt local_size, rest;

    rest = data.b0->size%mpi_size;
    local_size = data.b0->size/mpi_size + (rest > mpi_rank);
  
    ierr = VecCopy(state, rhs); CHKERRQ(ierr);
    ierr = VecScale(rhs, 0.5); CHKERRQ(ierr);
    ierr = VecGetArray(rhs, &rhs_array); CHKERRQ(ierr);
  
    for (PetscInt elem = 0; elem < local_size; ++elem) {
      OpOnBasisElems(env, data.b0->int_basis[elem], spin);
      rhs_array[elem] *= value;
    }
  
    ierr = VecRestoreArray(rhs, &rhs_array); CHKERRQ(ierr);

#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("SqzOpOnStateVector");
#endif
  
  return ierr;
}

/* --------------------------------------------------------------------------- */
// --------------------------------- Sqp Class ------------------------------- //
/* --------------------------------------------------------------------------- */

Sqp::Sqp(Environment& m_env,
	 Sq_data& m_data)
  : nspins{m_env.nspins},
    data{m_data},
    mpi_size{m_env.mpi_size},
    mpi_rank{m_env.mpi_rank},
    value{-1}
{}

/* --------------------------------------------------------------------------- */

void Sqp::OpOnBasisElems(Environment& env, PetscInt basis_elem, PetscInt spin) {

  PetscInt bitA, new_elem;
  value = -1;
  
  bitA = (basis_elem >> spin) & 0x1;
  if (!bitA) {
    new_elem = basis_elem ^ (1ull << spin);
    value = HamiltHelper::search_index(new_elem, data.b1->nspins);
  }
}

/* --------------------------------------------------------------------------- */

PetscErrorCode Sqp::OpOnStateVector(Environment& env, Vec& state, Vec& rhs, PetscInt spin) {

  PetscErrorCode ierr = 0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "SqpOpOnStateVector"};
#endif

    PetscScalar * state_array, * Jcoup;
    PetscInt * coup_elems;
    PetscInt ncol, istart, iend;

    ierr = VecSet(rhs,0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);

    ierr = VecGetOwnershipRange(state, &istart, &iend); CHKERRQ(ierr);
    
    ierr = PetscCalloc1(iend-istart, &coup_elems); CHKERRQ(ierr);
    ierr = PetscCalloc1(iend-istart, &Jcoup); CHKERRQ(ierr);
    ierr = VecGetArray(state, &state_array); CHKERRQ(ierr);

    ncol = 0;
    for (PetscInt elem = 0; elem < iend-istart; ++elem) {
      OpOnBasisElems(env, data.b0->int_basis[elem], spin);
      if (value != -1) {
    	coup_elems[ncol] = value;
    	Jcoup[ncol] = state_array[elem];
    	ncol += 1;
      }
    }
    ierr = VecSetValues(rhs,
			ncol,
			coup_elems,
			Jcoup,
			ADD_VALUES);
    CHKERRQ(ierr);
    
    ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);

    ierr = VecRestoreArray(state, &state_array); CHKERRQ(ierr);
    ierr = PetscFree(coup_elems); CHKERRQ(ierr);
    ierr = PetscFree(Jcoup); CHKERRQ(ierr);
    
#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("SqpOpOnStateVector");
#endif
  
  return ierr;
}
