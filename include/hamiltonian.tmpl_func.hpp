#ifndef __HAMILTONIAN_TMPL_FUNC_H
#define __HAMILTONIAN_TMPL_FUNC_H

/* --------------------------------------------------------------------------- */
// --------------------- Constructors/Destructors ---------------------------- //
/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Hamiltonian<L>::MatInit(Environment& env) {

  PetscErrorCode ierr = 0;
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "MatInit"};
    {
      Tools::ScopedTimer _timer_{env.tm, "Miscellaneous"};
#endif
      ierr = MatCreate(PETSC_COMM_WORLD, &hamilt); CHKERRQ(ierr);
      ierr = MatSetType(hamilt, MATMPIAIJ); CHKERRQ(ierr);
      ierr = MatSetSizes(hamilt, PETSC_DECIDE, PETSC_DECIDE, size, size); CHKERRQ(ierr);
      
#ifdef TIME_CODE
    }
    {
      Tools::ScopedTimer _timer_{env.tm, "prealloc"};
#endif
      ierr = prealloc(env); CHKERRQ(ierr);
#ifdef TIME_CODE
    }
    
    {
      Tools::ScopedTimer _timer_{env.tm, "MatZeroEntries"};
#endif
      ierr = MatZeroEntries(hamilt); CHKERRQ(ierr);
#ifdef TIME_CODE
    }
    
    {
      Tools::ScopedTimer _timer_{env.tm, "MatAssembly-init"};
#endif
      ierr = MatAssemblyBegin(hamilt, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(hamilt, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
#ifdef TIME_CODE
    }
    env.tm.PrintTimeInfoFunc("prealloc");
    env.tm.PrintTimeInfoFunc("Miscellaneous");
    env.tm.PrintTimeInfoFunc("MatZeroEntries");
    env.tm.PrintTimeInfoFunc("MatAssembly-init");
#endif
    
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("MatInit");
#endif
  return ierr;
}

/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Hamiltonian<L>::MatClean() {
  PetscErrorCode ierr = 0;
  ierr = MatDestroy(&hamilt); CHKERRQ(ierr);
  return ierr;
}

/* --------------------------------------------------------------------------- */

template<typename L>
Hamiltonian<L>::Hamiltonian (Environment& env,
			     Basis& m_basis,
			     L& m_lattice,
			     PetscReal m_J1,
			     PetscReal m_Delta1,
			     PetscReal m_J2,
			     PetscReal m_Delta2)
  : basis{m_basis},
    lattice{m_lattice},
    size{basis.size},
    J1{m_J1},
    Delta1{m_Delta1},
    J2{m_J2},
    Delta2{m_Delta2},
    mpi_size{env.mpi_size},
    mpi_rank{env.mpi_rank}
{
  PetscErrorCode ierr = 0;
  ierr = MatInit(env);
  if (ierr != 0)
    throw std::runtime_error("The matrix has not been properly initialized.\n");
}

/* --------------------------------------------------------------------------- */
// ------------------------- Building Functions ------------------------------ //
/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Hamiltonian<L>::build_disorder(Environment& env,
					      PetscReal* hi) {

  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend, local_index, tmp_elem1, tmp_elem2;
  PetscReal disorder_diag;

  ierr = MatGetOwnershipRange(hamilt, &Istart, &Iend); CHKERRQ(ierr);
  for (PetscInt elem = Istart; elem < Iend; ++elem) {
    local_index = elem - basis.global_start_index;

    tmp_elem1 = HamiltHelper::get_elem_diag(neigh_type::nn,basis.int_basis[local_index], basis.nspins, lattice);
    tmp_elem2 = HamiltHelper::get_elem_diag(neigh_type::nnn,basis.int_basis[local_index], basis.nspins, lattice);
    disorder_diag = Delta1*0.25 * tmp_elem1 + Delta2*0.25 * tmp_elem2 + 0.5 * HamiltHelper::get_disorder_elem(basis.int_basis[local_index], basis.nspins, hi);

    ierr = MatSetValue(hamilt,
		       elem,
		       elem,
		       disorder_diag,
		       INSERT_VALUES);
    CHKERRQ(ierr);  
  }
  
  ierr = MatAssemblyBegin(hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return (ierr);
}

/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Hamiltonian<L>::build_diag(Environment& env) {

  PetscErrorCode ierr = 0;
  PetscInt Istart, Iend, tmp_elem, local_index;
  PetscReal elem;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "build_diag"};
#endif
  
    ierr = MatGetOwnershipRange(hamilt, &Istart, &Iend); CHKERRQ(ierr);

    // Do not construct diagonal
    if (std::fabs(Delta1) < 1e-14 && std::fabs(Delta2) < 1e-14) {
      return(ierr);
    }

    // Construct nearest neighbours diagonal
    if (std::fabs(Delta1) >= 1e-14) {
      for (PetscInt i = Istart; i < Iend; ++i) {
	local_index = i - basis.global_start_index;
	  
#ifdef TIME_CODE
	{
	  Tools::ScopedTimer _timer_{env.tm, "get_coup_elems-d-nn"};
#endif
	  tmp_elem = HamiltHelper::get_elem_diag(neigh_type::nn,basis.int_basis[local_index], basis.nspins, lattice);
	    
	  if (tmp_elem == 0)
	    continue;
	    
	  elem = Delta1*0.25 * tmp_elem;
#ifdef TIME_CODE
	}
	{
	  Tools::ScopedTimer _timer_{env.tm, "MatSetValues-diag-nn"};
#endif
	  ierr = MatSetValue(hamilt,
			     i,
			     i,
			     elem,
			     ADD_VALUES);
	  CHKERRQ(ierr);
#ifdef TIME_CODE
	}
#endif
      }
    } // -- If (Delta1)

      // Construct next-nearest neighbours diagonal
    if (std::fabs(Delta2) >= 1e-14) {
      for (PetscInt i = Istart; i < Iend; ++i) {
	local_index = i - basis.global_start_index;
	  
#ifdef TIME_CODE
	{
	  Tools::ScopedTimer _timer_{env.tm, "get_coup_elems-d-nnn"};
#endif
	  tmp_elem = HamiltHelper::get_elem_diag(neigh_type::nnn,basis.int_basis[local_index], basis.nspins, lattice);

	  if (tmp_elem == 0)
	    continue;
	    
	  elem = Delta2*0.25 * tmp_elem;
#ifdef TIME_CODE
	}
	{
	  Tools::ScopedTimer _timer_{env.tm, "MatSetValues-diag-nnn"};
#endif
	  ierr = MatSetValue(hamilt,
			     i,
			     i,
			     elem,
			     ADD_VALUES);
	  CHKERRQ(ierr);
#ifdef TIME_CODE
	}
#endif
      }
    } // -- If (Delta2)      

#ifdef TIME_CODE
  }
  if (std::fabs(Delta1) >= 1e-14) {
    env.tm.PrintTimeInfoFunc("get_coup_elems-d-nn");
    env.tm.PrintTimeInfoFunc("MatSetValues-diag-nn");
  }
  if (std::fabs(Delta2) >= 1e-14) {
    env.tm.PrintTimeInfoFunc("get_coup_elems-d-nnn");
    env.tm.PrintTimeInfoFunc("MatSetValues-diag-nnn");
  }
  env.tm.PrintTimeInfoFunc("build_diag");
#endif
  return ierr;
}

/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Hamiltonian<L>::build_off_diag(Environment& env) {

  PetscErrorCode ierr = 0;
  PetscInt ncol, Istart, Iend, local_index;
  PetscScalar * Jcoup;
  PetscInt * coup_elems;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "build_off_diag"};
#endif
  // TODO: This may/should not be hardcoded
  PetscCalloc1(basis.nspins*3, &Jcoup);
  PetscCalloc1(basis.nspins*3, &coup_elems);
  
  // Construct nearest neighbours off-diagonal
  if (std::fabs(J1) >= 1e-14) {
    
    for (PetscInt i = 0; i < basis.nspins*3; ++i)
      Jcoup[i] = J1*0.5;
    
    ierr = MatGetOwnershipRange(hamilt, &Istart, &Iend); CHKERRQ(ierr);
    for (PetscInt elem = Istart; elem < Iend; ++elem) {
      local_index = elem - basis.global_start_index;
      
#ifdef TIME_CODE
      {
	Tools::ScopedTimer _timer_{env.tm, "get_coup_elems-o-nn"};
#endif
	HamiltHelper::get_coup_elems(neigh_type::nn, basis, local_index, lattice, ncol, coup_elems);
#ifdef TIME_CODE
      }
      {
	Tools::ScopedTimer _timer_{env.tm, "MatSetValues-off-nn"};
#endif
	ierr = MatSetValues(hamilt,
			    1, /* number of rows */
			    &elem, /* row(s) */
			    ncol, /* number of columns */
			    coup_elems, /* columns */
			    Jcoup, /* values to fill in with */
			    ADD_VALUES);
	CHKERRQ(ierr);
#ifdef TIME_CODE
      }
#endif
    }
  } // -- if J1

  // Construct next-nearest neighbours off-diagonal
  if (std::fabs(J2) >= 1e-14) {
   
    for (PetscInt i = 0; i < basis.nspins*3; ++i)
      Jcoup[i] = J2*0.5;
    
    ierr = MatGetOwnershipRange(hamilt, &Istart, &Iend); CHKERRQ(ierr);
    for (PetscInt elem = Istart; elem < Iend; ++elem) {
      local_index = elem - basis.global_start_index;
      
#ifdef TIME_CODE
      {
	Tools::ScopedTimer _timer_{env.tm, "get_coup_elems-o-nnn"};
#endif
	HamiltHelper::get_coup_elems(neigh_type::nnn, basis, local_index, lattice, ncol, coup_elems);
#ifdef TIME_CODE
      }
      {
	Tools::ScopedTimer _timer_{env.tm, "MatSetValues-off-nnn"};
#endif
	ierr = MatSetValues(hamilt,
			    1, /* number of rows */
			    &elem, /* row(s) */
			    ncol, /* number of columns */
			    coup_elems, /* columns */
			    Jcoup, /* values to fill in with */
			    ADD_VALUES);
	CHKERRQ(ierr);
#ifdef TIME_CODE
      }
#endif
    }
  } // -- if J2

  ierr = PetscFree(Jcoup);
  ierr = PetscFree(coup_elems);
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "MatAssembly-offdiag"};
#endif

    if (env.disorder_flg == false) {
      ierr = MatAssemblyBegin(hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    } else {
      ierr = MatAssemblyBegin(hamilt,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(hamilt,MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr);    
    }
    
#ifdef TIME_CODE
  }
  }
  if (std::fabs(Delta1) >= 1e-14) {
    env.tm.PrintTimeInfoFunc("get_coup_elems-o-nn");
    env.tm.PrintTimeInfoFunc("MatSetValues-off-nn");
  }
  if (std::fabs(Delta2) >= 1e-14) {
    env.tm.PrintTimeInfoFunc("get_coup_elems-o-nnn");
    env.tm.PrintTimeInfoFunc("MatSetValues-off-nnn");
  }
  env.tm.PrintTimeInfoFunc("MatAssembly-offdiag");
  env.tm.PrintTimeInfoFunc("build_off_diag");
#endif
  return ierr;
}

/* --------------------------------------------------------------------------- */
// ---------------------------- Extra Functions ------------------------------ //
/* --------------------------------------------------------------------------- */

template<typename L>
void HamiltHelper::get_coup_elems(neigh_type nt,
				  Basis& basis,
				  PetscInt elem,
				  L& lattice,
				  PetscInt& ncol,
				  PetscInt* coup_elems) {


  switch(nt) {

    // Go over all nearest neighbours of spin i
  case neigh_type::nn:
    ncol = 0;
    for (PetscInt i = 0; i < basis.nspins; ++i) {
      for (PetscInt neigh = 1; neigh <= lattice.num_nnXsite(i); ++neigh) {
	
	PetscInt offset = (lattice.num_nn() + 1) * i;
	short int spin = lattice.nn[(offset+neigh)];
	get_coup_elems_AUX(basis, elem, i, spin, ncol, coup_elems);
      }
    }
    break;

    // Go over all next-nearest neighbours of spin i
  case neigh_type::nnn:
    ncol = 0;
    for (PetscInt i = 0; i < basis.nspins; ++i) {
      for (PetscInt neigh = 1; neigh <= lattice.num_nnnXsite(i); ++neigh) {
	
  	PetscInt offset = (lattice.num_nnn() + 1) * i;
  	short int spin = lattice.nnn[(offset+neigh)];
  	get_coup_elems_AUX(basis, elem, i, spin, ncol, coup_elems);   
      }
    }
    break;
  }
}

/* --------------------------------------------------------------------------- */

template<typename L>
PetscInt HamiltHelper::get_elem_diag(neigh_type nt,
				     PetscInt basis_elem,
				     PetscInt nspins,
				     L& lattice) {

  PetscInt diag_elem = 0;

  switch (nt) {

    // Go over all nearest neighbours of spin i
  case neigh_type::nn:
    for (PetscInt i = 0; i < nspins; ++i) {
      for (PetscInt neigh = 1; neigh <= lattice.num_nnXsite(i); ++neigh) {
	
	PetscInt offset = (lattice.num_nn() + 1) * i;
	short int spin = lattice.nn[(offset+neigh)];
	get_elem_diag_AUX(basis_elem, i, spin, diag_elem);
      } 
    }
    break;
    
    // Go over all nearest neighbours of spin i
  case neigh_type::nnn:
    for (PetscInt i = 0; i < nspins; ++i) {
      for (PetscInt neigh = 1; neigh <= lattice.num_nnnXsite(i); ++neigh) {
	
	PetscInt offset = (lattice.num_nnn() + 1) * i;
	short int spin = lattice.nnn[(offset+neigh)];
	get_elem_diag_AUX(basis_elem, i, spin, diag_elem);
      } 
    }
    break;
  }
  return diag_elem;
}

/* --------------------------------------------------------------------------- */
// ---------------------------- Prealloc Functions --------------------------- //
/* --------------------------------------------------------------------------- */

template<typename L>
PetscErrorCode Hamiltonian<L>::prealloc(Environment& env) {

  PetscErrorCode ierr = 0;
  PetscInt local_size = PETSC_DECIDE;
  PetscInt rest = size%mpi_size;
  PetscInt ncol, elem_diag = 0;
  PetscInt * d_nnz, * o_nnz, * coup_elems;
  PetscInt global_start_index, global_end_index, local_index;

  ierr = PetscSplitOwnership(MPI_COMM_WORLD, &local_size, &size); CHKERRQ(ierr);
  global_start_index = local_size*mpi_rank + (rest <= mpi_rank)*rest;
  global_end_index = local_size*(mpi_rank+1) + (rest <= mpi_rank)*rest;

  PetscCalloc1(local_size, &d_nnz);
  PetscCalloc1(local_size, &o_nnz);

  ncol = 0;

  // TODO: This may/should not be hardcoded
  PetscCalloc1(basis.nspins*3, &coup_elems);

  if (std::fabs(J1) >= 1e-14) {
    for (PetscInt elem = global_start_index; elem < global_end_index; ++elem) {
      
      local_index = elem - basis.global_start_index;

#ifdef TIME_CODE
      {
	Tools::ScopedTimer _timer_{env.tm, "prealloc-J1"};
#endif
	  
	HamiltHelper::get_coup_elems(neigh_type::nn, basis, local_index, lattice, ncol, coup_elems);
	HamiltHelper::prealloc_info(Delta1,
				    elem,
				    global_start_index,
				    global_end_index,
				    coup_elems,
				    ncol,
				    elem_diag,
				    d_nnz,
				    o_nnz);
	
#ifdef TIME_CODE
      }
#endif
    } // -- for
  } // -- if J1
  
  if (std::fabs(J2) >= 1e-14) {
    ncol = 0; // VERY IMPORTANT: It's the way to avoid preallocating for the J2 part if not needed
    for (PetscInt elem = global_start_index; elem < global_end_index; ++elem) {
      local_index = elem - basis.global_start_index;
      
#ifdef TIME_CODE
      {
	Tools::ScopedTimer _timer_{env.tm, "prealloc-J2"};
#endif
	  
	HamiltHelper::get_coup_elems(neigh_type::nnn, basis, local_index, lattice, ncol, coup_elems);
	HamiltHelper::prealloc_info(Delta2,
				    elem,
				    global_start_index,
				    global_end_index,
				    coup_elems,
				    ncol,
				    elem_diag,
				    d_nnz,
				    o_nnz);
#ifdef TIME_CODE
      }
#endif
    } // -- for
  } // -- if J2
  
  ierr = PetscFree(coup_elems);

  //#ifndef DISORDER
  // Add one for the diagonal element if at least one Delta is not zero
  //if ((std::fabs(Delta1) >= 1e-14) | (std::fabs(Delta2) >= 1e-14))
    //#endif
    for (PetscInt i = 0; i < local_size; ++i)
      d_nnz[i] += 1;
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "PetscPreallocation"};
#endif
    
    ierr = MatMPIAIJSetPreallocation(hamilt, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);
#ifdef TIME_CODE
  }
  if (std::fabs(J1) >= 1e-14) env.tm.PrintTimeInfoFunc("prealloc-J1");
  if (std::fabs(J2) >= 1e-14) env.tm.PrintTimeInfoFunc("prealloc-J2");
  env.tm.PrintTimeInfoFunc("PetscPreallocation");
#endif
  ierr = PetscFree(d_nnz);
  ierr = PetscFree(o_nnz);
  return ierr;
}

#endif
