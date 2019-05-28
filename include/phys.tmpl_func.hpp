/* --------------------------------------------------------------------------- */
// ------------------------ Dynamical Structure Factor ----------------------- //
/* --------------------------------------------------------------------------- */

template<typename L, template<typename> class SQ>
PetscScalar Phys::DynStructFactor(Environment& env,
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
    PetscReal dsf2;
    std::array<PetscReal,2> q;
  
    Vec rhs, eig_vec;
    ierr = MatCreateVecs(data.h1->hamilt,&eig_vec,NULL); CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, data.b1->size, &rhs); CHKERRQ(ierr);
    // ierr = VecDuplicate(state, &rhs); CHKERRQ(ierr);

    EPS DSF_solver;
    PetscInt nconv;
    PetscScalar eig_val;
    PetscReal error;
    std::vector<PetscScalar> c0, eigenvals;
    PetscInt offset;

#ifdef TIME_CODE
    {
      Tools::ScopedTimer _timer_{env.tm, "SolverDSF"};
#endif
      
    // Initiate solver and solve
    ierr = Solver::SolverInit(DSF_solver, data.h1->hamilt, data.nev, data.ncv, data.mpd); CHKERRQ(ierr);
    ierr = Solver::solve_lanczos(env, DSF_solver, nconv); CHKERRQ(ierr);
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

    dsf_fp = fopen(dsf_filename.str().c_str(), "a");
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

    for (PetscInt iter_q = 0; iter_q < data.nQ; ++ iter_q) {

      q = data.lat->get_q((*data.qi)[iter_q]);

      ierr = PetscFPrintf(PETSC_COMM_WORLD, dsf_fp, "# iter_q: %d\n", (*data.qi)[iter_q]); CHKERRQ(ierr);
    
      // Build Sq operator
      Sq_data<L> sq_data{data.b0, data.b1, data.lat, q};
      SQ<L> sq_oper{env, sq_data};
    
      // Compute initial vector
      ierr = sq_oper.OpOnStateVector(env, state, rhs); CHKERRQ(ierr);

#ifdef TIME_CODE
      {
	Tools::ScopedTimer _timer_{env.tm,"CoeffDSF"};
#endif
      // Get eigenvectors and do the product with the rhs vector
      for (PetscInt i = 0; i < data.nev; ++i) {
	ierr = Solver::solution(DSF_solver, i, eig_val, eig_vec, error); CHKERRQ(ierr);
	eigenvals.push_back(eig_val);
	ierr = VecDot(rhs, eig_vec, &dsf); CHKERRQ(ierr); // --> VecDot(x,y) = y^H x
	dsf2 = std::norm<PetscReal>(dsf);
	c0.push_back(dsf2);
      }
#ifdef TIME_CODE
      }
#endif

      // Write to file
      offset = data.nev * iter_q;
      for (PetscInt i = offset; i < offset + data.nev; ++i) {
	ierr = PetscFPrintf(PETSC_COMM_WORLD,
			    dsf_fp,
			    "%.13e  %.13e\n",
			    PetscRealPart(eigenvals[i] - gs),
			    PetscRealPart(c0[i])); CHKERRQ(ierr);
      } // -- for nconv
    } // -- for iter_q

    ierr = Solver::SolverClean(DSF_solver); CHKERRQ(ierr);
    ierr = VecDestroy(&rhs); CHKERRQ(ierr);
    ierr = VecDestroy(&eig_vec); CHKERRQ(ierr);
    fclose(dsf_fp);

#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("CoeffDSF");
  env.tm.PrintTimeInfoFunc("SolverDSF");
  env.tm.PrintTimeInfoFunc("DynStructFactor");
#endif

  return dsf;
}
