/* --------------------------------------------------------------------------- */
// ------------------------ Magnetic Order Parameter ------------------------- //
/* --------------------------------------------------------------------------- */

template<typename L>
PetscScalar Phys::SquaredSubLattMagAF(Environment& env,
				      Basis& basis,
				      AF<L>& sublat,
				      Vec& state) {

  PetscScalar mAF = 0.0;
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "MagAF"};
#endif
    
    SiSj cij{env,basis};
    
    PetscScalar exp_val = 0.0;
    PetscInt spin_i, spin_j;
    
    for (PetscInt lat = 0; lat < sublat.get_size(); ++lat) {
      for (PetscInt i = 0; i < PetscInt(sublat.get_sl(lat).size()); ++i) {
	spin_i = sublat.get_sl(lat)[i];
	for (PetscInt j = 0; j < PetscInt(sublat.get_sl(lat).size()); ++j) {
	  spin_j = sublat.get_sl(lat)[j];
	  exp_val += cij.ExpectVal(env,spin_i, spin_j, state);
	}
      }
    }
    
    mAF = exp_val * 8.0 / PetscScalar(env.nspins) / (PetscScalar(env.nspins + 4.0));
    
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("MagAF");
#endif
  
  return mAF;
}

/* --------------------------------------------------------------------------- */

template<typename L>
PetscScalar Phys::SquaredSubLattMagSTR(Environment& env,
				       Basis& basis,
				       Striped<L>& sublat,
				       Vec& state) {

  PetscScalar mSTR = 0.0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "MagSTR"};
#endif
    
    SiSj cij{env,basis};

    PetscScalar exp_val = 0.0;
    PetscInt spin_i, spin_j;

    for (PetscInt str = 0; str < 2; ++str)
      for (PetscInt lat = 0; lat < 2; ++lat)
	for (PetscInt i = 0; i < PetscInt(sublat.get_sl(str*2+lat).size()); ++i) {
	  spin_i = sublat.get_sl(str*2+lat)[i];
	  for (PetscInt j = 0; j < PetscInt(sublat.get_sl(str*2+lat).size()); ++j) {
	    spin_j = sublat.get_sl(str*2+lat)[j];
	    exp_val += cij.ExpectVal(env,spin_i, spin_j, state);
	  }
	}

    mSTR = exp_val * 4.0 / PetscScalar(env.nspins) / (PetscScalar(env.nspins + 4.0));

#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("MagSTR");
#endif
  
  return mSTR;
}

/* --------------------------------------------------------------------------- */
// ------------------------ Dynamical Structure Factor ----------------------- //
/* --------------------------------------------------------------------------- */

template<typename L, typename SQ>
PetscErrorCode Phys::DynStructFactor(Environment& env,
				     Phys::DSF_data<L>& data,
				     PetscScalar gs,
				     Vec& state,
				     PetscInt dis_iter) {

  PetscScalar dsf;
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "DynStructFactor"};
#endif
    
    PetscErrorCode ierr = 0;
    
    Vec rhs, eig_vec;
    ierr = MatCreateVecs(data.h1->hamilt,&eig_vec,NULL); CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, data.b1->size, &rhs); CHKERRQ(ierr);
    //ierr = VecDuplicate(state, &rhs); CHKERRQ(ierr);

    EPS DSF_solver;
    PetscInt nconv;
    PetscScalar eig_val;
    PetscReal error;

#ifdef TIME_CODE
    {
      Tools::ScopedTimer _timer_{env.tm, "SolverDSF"};
#endif
      
      // Initiate solver and solve
      ierr = Solver::SolverInit(DSF_solver, data.h1->hamilt, data.nev, data.ncv, data.mpd); CHKERRQ(ierr);
      ierr = Solver::solve(env, DSF_solver, nconv); CHKERRQ(ierr);
#ifdef TIME_CODE
    }
#endif
  
    // ------------------------------- FILE -------------------------------------
    std::ostringstream dsf_filename;
    FILE * dsf_fp;

    dsf_filename << data.path
		 << "/dsf"
		 << "_d" << dis_iter
		 << ".dat";

    ierr = PetscFOpen(PETSC_COMM_SELF,dsf_filename.str().c_str(), "w", &dsf_fp); CHKERRQ(ierr);
    if (!dsf_fp) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Cannot open file");

    ierr = PetscFPrintf(PETSC_COMM_WORLD,
			dsf_fp, "#------------------------------- Parameters -------------------------------------\n"); CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,
			dsf_fp, "#%7s:   %d  \n", " nspins", env.nspins); CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,
			dsf_fp, "#%7s:   %d  \n", " nev", data.nev); CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,
			dsf_fp, "#%7s:   %.16f  \n", " GS", PetscRealPart(gs)); CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,
			dsf_fp, "#%7s:   %.5f  %4s:   %.5f  %4s:   %.5f  %4s:   %.5f\n",
			"J1", data.h0->hamilt_J1(), "D1", data.h0->hamilt_D1(), "J2", data.h0->hamilt_J2(), "D2", data.h0->hamilt_D2()); CHKERRQ(ierr);

    if (data.disorder == PETSC_TRUE) {
      for (PetscInt i = 0; i < env.nspins; ++i)
	ierr = PetscFPrintf(PETSC_COMM_WORLD,
			    dsf_fp, "# hi[%4d]: %.16f\n", i, data.hi[i]); CHKERRQ(ierr);
    }
  
    ierr = PetscFPrintf(PETSC_COMM_WORLD,
			dsf_fp, "#--------------------------------------------------------------------------------\n#\n#\n"); CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,
			dsf_fp, "#----- Energy -----  ------ Factor -----\n"); CHKERRQ(ierr);
    // ------------------------------- FILE -------------------------------------

    // Build Sq operator
    Sq_data sq_data{data.b0, data.b1};
    SQ sq_oper{env, sq_data};
    std::vector<PetscScalar> c0, eigenvals;
    
    for (PetscInt iter_s = 0; iter_s < env.nspins; ++iter_s) {
        
      // Compute initial vector
      ierr = sq_oper.OpOnStateVector(env, state, rhs, iter_s); CHKERRQ(ierr);

#ifdef TIME_CODE
      {
	Tools::ScopedTimer _timer_{env.tm,"CoeffDSF"};
#endif
	// Get eigenvectors and do the product with the rhs vector
	for (PetscInt i = 0; i < data.nev; ++i) {
	  ierr = Solver::solution(DSF_solver, i, eig_val, eig_vec, error); CHKERRQ(ierr);
	  eigenvals.push_back(eig_val);
	  ierr = VecDot(rhs, eig_vec, &dsf); CHKERRQ(ierr); // --> VecDot(x,y) = y^H x
	  c0.push_back(dsf);
	}
#ifdef TIME_CODE
      }
#endif
      
    } // -- for iter_s

    // Write to file
    PetscInt offset;
    for (PetscInt spin = 0; spin < env.nspins; ++spin){
      ierr = PetscFPrintf(PETSC_COMM_WORLD, dsf_fp, "# spin: %d\n", spin); CHKERRQ(ierr);
      offset = spin*data.nev;
      for (PetscInt i = 0; i < data.nev; ++i) {
	ierr = PetscFPrintf(PETSC_COMM_WORLD,
			    dsf_fp,
			    "%.13e  %.13e\n",
			    PetscRealPart(eigenvals[i+offset] - gs),
			    PetscRealPart(c0[i+offset])); CHKERRQ(ierr);
      }
    }
    
    ierr = Solver::SolverClean(DSF_solver); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs); CHKERRQ(ierr);
    ierr = VecDestroy(&eig_vec); CHKERRQ(ierr);
    ierr = PetscFClose(PETSC_COMM_SELF, dsf_fp); CHKERRQ(ierr);

#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("CoeffDSF");
  env.tm.PrintTimeInfoFunc("SolverDSF");
  env.tm.PrintTimeInfoFunc("DynStructFactor");
#endif

  return ierr;
}
