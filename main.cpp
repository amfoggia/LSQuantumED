#include "sublattice.hpp"
#include "time_monitor.hpp"
#include "random_disorder.hpp"
#include "parse_args.hpp"
#include "print_args.hpp"
#include "phys.hpp"
#include <cstdlib>
#include <petsctime.h>

static char help[] = "1/2-spin system governed by antiferromagnetic Heisenberg Hamiltonian with nn and nnn interactions, with constant total magnetization\n";

int main(int argc, char * argv[]) {

  EPS solver;
  
  PetscErrorCode ierr = 0;
  PetscLogStage stage1, stage2, stage3, stage4a, stage4b, stage7;
  PetscInt nconv, spins;
  PetscScalar real_eig;
  Vec real_vec;
  PetscReal error;
  PetscScalar exp_val, magAF = 0.0;

  if (argc < 2) {
    std::cout << "At leats one argument is needed. "
  	      << "Number of total spins need to be provided after the executable's name.\n";
    return 1;
  }
  
  std::istringstream iss(argv[1]);
  if ( !(iss >> spins) ) {
    std::cout << "Number of total spins need to be provided after the executable's name.\n";
    return 1;
  }
  
  Environment env{argc,argv,spins,help};
  
  parse_args();
  print_args();

#ifdef DISORDER
  Tools::RandomDisorder rd{spins,reprod,min_dis,max_dis,rep};
#endif
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "1- Creating basis\n");
  
  PetscLogStageRegister("Create Basis", &stage1);
  PetscLogStagePush(stage1);
  Basis b{env,0};
  PetscLogStagePop();
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "2- Creating lattice\n"); CHKERRQ(ierr);

  PetscLogStageRegister("Create Lattice", &stage2);
  PetscLogStagePush(stage2);
  chain1D lattice{env};
  //  AF<2,square2D> sub_af{lattice};
  PetscLogStagePop();
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "3- Creating Hamiltonian\n"); CHKERRQ(ierr);

  PetscLogStageRegister("Create Hamilt", &stage3);
  PetscLogStagePush(stage3);
  Hamiltonian<chain1D> h{env,b,lattice,J1,Delta1,J2,Delta2};
  PetscLogStagePop();

  ierr = PetscPrintf(MPI_COMM_WORLD, "4a- Building diagonal part\n"); CHKERRQ(ierr);

  PetscLogStageRegister("Diag", &stage4a);
  PetscLogStagePush(stage3);
  ierr = h.build_diag(env); CHKERRQ(ierr);
  PetscLogStagePop();
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "4b- Building offdiagonal part\n"); CHKERRQ(ierr);

  PetscLogStageRegister("Offdiag", &stage4b);
  PetscLogStagePush(stage4b);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);
  PetscLogStagePop();

#ifndef DISORDER
  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
#endif
    
  //  ierr = PetscPrintf(MPI_COMM_WORLD, "4c- Building disorder part\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     " iter          k          ||Ax-kx||/||kx||\n"
		     " ----    -------------   ------------------\n"); CHKERRQ(ierr);

  // DSF input data
  std::vector<PetscInt> tmp;
  for (PetscInt i = 0; i < spins; ++i)
    tmp.push_back(i);

  PetscScalar dynstructfactor;
  Phys::DSF_data<chain1D> dsf_data{};
  dsf_data.b0 = &b;
  dsf_data.b1 = &b;
  dsf_data.lat = &lattice;
  dsf_data.nQ = spins;
  dsf_data.qi = &tmp;
  dsf_data.nev = 10;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
#ifdef DISORDER
  dsf_data.disorder = PETSC_TRUE;
#else
  dsf_data.disorder = PETSC_FALSE;
#endif
  dsf_data.path = "./";

#ifndef DISORDER
  ierr = MatCreateVecs(h.hamilt,&real_vec,NULL); CHKERRQ(ierr);
#endif

#ifdef DISORDER
  for (int repet = 0; repet < rep; ++repet) {
    rd.hi_init(repet);
    ierr = h.build_disorder(env,rd.hi); CHKERRQ(ierr);
    dsf_data.hi = rd.hi;
#else
    int repet = 0;
    dsf_data.hi = NULL;
#endif
    dsf_data.h0 = &h;
    dsf_data.h1 = &h;
    ierr = Solver::SolverInit(solver,h.hamilt,1); CHKERRQ(ierr);
    ierr = Solver::solve_lanczos(env,solver,nconv); CHKERRQ(ierr);
#ifdef DISORDER
    ierr = MatCreateVecs(h.hamilt,&real_vec,NULL); CHKERRQ(ierr);
#endif
    if (nconv > 0) {      
      for (PetscInt i = 0; i < 1; ++i) {
	TIME_IT(env.tm, "Get Eigen", ierr = Solver::solution(solver,i,real_eig,real_vec,error)); CHKERRQ(ierr);
	ierr = Solver::solution(solver,i,real_eig,real_vec,error); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"  %d     %12f       %12g\n",repet,PetscRealPart(real_eig),(double)error);
	CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);
    } // -- if nconv

#ifdef TIME_CODE
    env.tm.PrintTimeInfoFunc("Get Eigen");
#endif

    // exp_val = 0.0; 
    // exp_val = Phys::SquaredSubLattMagAF(env, b, std::array<std::vector<PetscInt>,2>{sub_af.get_sl(0),sub_af.get_sl(1)}, real_vec);
    // magAF += exp_val;
    
    // exp_val = 0.0;
    // Striped<4,square2D> sub_str{lattice};
    // exp_val = Phys::SquaredSubLattMagSTR(env, b, std::array<std::array<std::vector<PetscInt>,2>,2>{sub_str.get_sl(0),sub_str.get_sl(1),sub_str.get_sl(2),sub_str.get_sl(3)}, real_vec);
    // ierr = PetscPrintf(PETSC_COMM_WORLD, "exp_val: %.10f\n", exp_val);

    // PetscPrintf(PETSC_COMM_WORLD, "magAF: %.10f\n", magAF/PetscReal(rep)); CHKERRQ(ierr);


    ierr = PetscPrintf(MPI_COMM_WORLD, "8- Phys quantities\n"); CHKERRQ(ierr);

#ifndef DISORDER
    PetscLogStageRegister("Phys quantities", &stage7);
    PetscLogStagePush(stage7);
#endif
    dynstructfactor = Phys::DynStructFactor<chain1D,Sqz>(env, dsf_data, real_eig, real_vec, repet);
#ifndef DISORDER
    PetscLogStagePop();
#endif

    ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
#ifdef DISORDER
    ierr = VecDestroy(&real_vec); CHKERRQ(ierr);
  } // -- for repet
#endif

#ifndef DISORDER
  ierr = VecDestroy(&real_vec); CHKERRQ(ierr);
#endif
  ierr = PetscPrintf(MPI_COMM_WORLD, "9- Cleaning\n"); CHKERRQ(ierr);

#ifdef TIME_CODE
  env.tm.PrintTimeInfoFull();
#endif

  return ierr;
}
