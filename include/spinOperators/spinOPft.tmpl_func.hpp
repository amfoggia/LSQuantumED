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
{
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{m_env.tm, "SqzConstructor"};
#endif
    
    exp = compute_exp();
    exp.shrink_to_fit();

#ifdef TIME_CODE
  }
  //m_env.tm.PrintTimeInfoFunc("SqzConstructor");
#endif
  
}

/* --------------------------------------------------------------------------- */

template<typename L>
std::vector<PetscScalar>& Sqz<L>::compute_exp() {
  
  PetscReal qxr;
  PetscComplex power;
  for (PetscInt i = 0; i < nspins; ++i) {
    qxr = data.qi[0]*PetscReal(data.lat->get_r(i)[0]);
    qxr += data.qi[1]*PetscReal(data.lat->get_r(i)[1]);
    power = PETSC_i * qxr;
    exp.push_back(std::exp(power));
  }
  return exp;
}

/* --------------------------------------------------------------------------- */

template<typename L>
void Sqz<L>::OpOnBasisElems(Environment& env, PetscInt basis_elem) {

  PetscInt bitA;
  value = 0.0;

  for(PetscInt i = 0; i < nspins; ++i) {
    bitA = (basis_elem >> i) & 0x1;
    if (bitA)
      value = value + exp[i];
    else
      value = value - exp[i];
  }
}

/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Sqz<L>::OpOnStateVector(Environment& env, Vec& state, Vec& rhs) {

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
    ierr = VecScale(rhs, 0.5/std::sqrt(nspins)); CHKERRQ(ierr);
    ierr = VecGetArray(rhs, &rhs_array); CHKERRQ(ierr);
  
    for (PetscInt elem = 0; elem < local_size; ++elem) {
      OpOnBasisElems(env, data.b0->int_basis[elem]);
      rhs_array[elem] *= value;
    }
  
    ierr = VecRestoreArray(rhs, &rhs_array); CHKERRQ(ierr);

#ifdef TIME_CODE
  }
  //env.tm.PrintTimeInfoFunc("SqzOpOnStateVector");
#endif
  
  return ierr;
}
