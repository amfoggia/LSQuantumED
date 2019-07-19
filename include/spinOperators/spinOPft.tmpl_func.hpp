/* --------------------------------------------------------------------------- */
// --------------------------------- Sqz Class ------------------------------- //
/* --------------------------------------------------------------------------- */

template<typename L>
Sqz<L>::Sqz(Environment& m_env,
	    Sq_data<L>& m_data)
  : nspins{m_env.nspins},
    data{m_data},
    mpi_size{m_env.mpi_size},
    mpi_rank{m_env.mpi_rank},
    value{0.0}
{}

/* --------------------------------------------------------------------------- */

template<typename L>
void Sqz<L>::OpOnBasisElems(Environment& env, PetscInt basis_elem, PetscInt spin) {

  PetscInt bitA;
  value = 1.0;
  
  bitA = (basis_elem >> spin) & 0x1;
  if (!bitA)
    value *= -1.0;
}

/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Sqz<L>::OpOnStateVector(Environment& env, Vec& state, Vec& rhs, PetscInt spin) {

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
