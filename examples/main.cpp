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
  PetscLogStage stage1, stage2, stage3, stage4a, stage4b, stage4c, stage5, stage5a, stage5b, stage6, stage7;
  PetscInt nconv, spins;
  PetscScalar real_eig;
  Vec real_vec;
  PetscReal error;
  PetscScalar exp_val, magAF = 0.0;
  
  Environment env{argc,argv,help};
  
  parse_args();
  print_args();

  PetscLogStageRegister("Create Basis", &stage1);
  PetscLogStageRegister("Create Lattice", &stage2);
  PetscLogStageRegister("Create Hamilt", &stage3);
  PetscLogStageRegister("Diag", &stage4a);
  PetscLogStageRegister("Offdiag", &stage4b);
  PetscLogStageRegister("FinaleAssembly", &stage4c);
  PetscLogStageRegister("SolvingProblem", &stage5);
  PetscLogStageRegister("Build disorder", &stage5a);
  PetscLogStageRegister("Solve disorder", &stage5b);
  PetscLogStageRegister("Get GS", &stage6);
  PetscLogStageRegister("Phys quantities", &stage7);

  Tools::RandomDisorder rd{env.nspins,reprod,min_dis,max_dis,rep};

  ierr = PetscPrintf(MPI_COMM_WORLD, "1- Creating basis\n");  
  PetscLogStagePush(stage1);
  Basis b{env,0};
  PetscLogStagePop();
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "2- Creating lattice\n"); CHKERRQ(ierr);
  PetscLogStagePush(stage2);
  chain1D lattice{env};
  //square2D lattice{env,4,4};
  //AF<2,square2D> sub_af{lattice};
  PetscLogStagePop();
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "3- Creating Hamiltonian\n"); CHKERRQ(ierr);
  PetscLogStagePush(stage3);
  Hamiltonian<chain1D> h{env,b,lattice,J1,Delta1,J2,Delta2};
  //Hamiltonian<square2D> h{env,b,lattice,J1,Delta1,J2,Delta2};
  PetscLogStagePop();

  ierr = PetscPrintf(MPI_COMM_WORLD, "4a- Building diagonal part\n"); CHKERRQ(ierr);
  PetscLogStagePush(stage4a);
  ierr = h.build_diag(env); CHKERRQ(ierr);
  PetscLogStagePop();
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "4b- Building offdiagonal part\n"); CHKERRQ(ierr);
  PetscLogStagePush(stage4b);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);
  PetscLogStagePop();

  if (env.disorder_flg == false) {
    ierr = PetscPrintf(MPI_COMM_WORLD, "4c- Final Assembly\n"); CHKERRQ(ierr);
    PetscLogStagePush(stage4c);
    ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    PetscLogStagePop();
  }
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "5- Solving Problem\n"); CHKERRQ(ierr);
  PetscLogStagePush(stage5);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     " iter          k          ||Ax-kx||/||kx||\n"
		     " ----    -------------   ------------------\n"); CHKERRQ(ierr);
  
  PetscScalar dynstructfactor;
  Phys::DSF_data<chain1D> dsf_data{};
  //Phys::DSF_data<square2D> dsf_data{};
  dsf_data.b0 = &b;
  dsf_data.b1 = &b;
  dsf_data.lat = &lattice;
  dsf_data.nev = 10;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
  dsf_data.disorder = env.disorder_flg;
  dsf_data.path = "./";

  if (env.disorder_flg == false)
    ierr = MatCreateVecs(h.hamilt,&real_vec,NULL); CHKERRQ(ierr);

  for (int repet = 0; repet < rep; ++repet) {
    if (env.disorder_flg == true) {
      rd.hi_init(repet);
      
      PetscLogStagePush(stage5a);
      ierr = h.build_disorder(env,rd.hi); CHKERRQ(ierr);
      PetscLogStagePop();
      
      dsf_data.hi = rd.hi;
    }
    if (env.disorder_flg == false)
      dsf_data.hi = NULL;
    dsf_data.h0 = &h;
    dsf_data.h1 = &h;
    
    PetscLogStagePush(stage5b);
    ierr = Solver::SolverInit(solver,h.hamilt,1); CHKERRQ(ierr);
    ierr = Solver::solve_lanczos(env,solver,nconv); CHKERRQ(ierr);
    PetscLogStagePop();
    
    if (env.disorder_flg == true)
      ierr = MatCreateVecs(h.hamilt,&real_vec,NULL); CHKERRQ(ierr);
    
    PetscLogStagePush(stage6);
    if (nconv > 0) {      
      for (PetscInt i = 0; i < 1; ++i) {
	TIME_IT(env.tm, "Get Eigen", ierr = Solver::solution(solver,i,real_eig,real_vec,error)); CHKERRQ(ierr);
	ierr = Solver::solution(solver,i,real_eig,real_vec,error); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"  %d     %12f       %12g\n",repet,PetscRealPart(real_eig),(double)error);
	CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);
    } // -- if nconv
    PetscLogStagePop(); // --> stage6 - GS

#ifdef TIME_CODE
    env.tm.PrintTimeInfoFunc("Get Eigen");
#endif
    
    ierr = PetscPrintf(MPI_COMM_WORLD, "8- Phys quantities\n"); CHKERRQ(ierr);
    // exp_val = 0.0; 
    // exp_val = Phys::SquaredSubLattMagAF(env, b, std::array<std::vector<PetscInt>,2>{sub_af.get_sl(0),sub_af.get_sl(1)}, real_vec);
    // magAF += exp_val;
    
    // exp_val = 0.0;
    // Striped<4,square2D> sub_str{lattice};
    // exp_val = Phys::SquaredSubLattMagSTR(env, b, std::array<std::array<std::vector<PetscInt>,2>,2>{sub_str.get_sl(0),sub_str.get_sl(1),sub_str.get_sl(2),sub_str.get_sl(3)}, real_vec);
    // ierr = PetscPrintf(PETSC_COMM_WORLD, "exp_val: %.10f\n", exp_val);
    
    // PetscPrintf(PETSC_COMM_WORLD, "magAF: %.10f\n", magAF/PetscReal(rep)); CHKERRQ(ierr);
    
    PetscLogStagePush(stage7);
    dynstructfactor = Phys::DynStructFactor<chain1D,Sqz>(env, dsf_data, real_eig, real_vec, repet);
    //dynstructfactor = Phys::DynStructFactor<square2D,Sqz>(env, dsf_data, real_eig, real_vec, repet);
    PetscLogStagePop(); // --> stage7 - Phys quantities
    
    ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
    if (env.disorder_flg == true)
      ierr = VecDestroy(&real_vec); CHKERRQ(ierr);
  } // -- for repet
  PetscLogStagePop(); // --> stage5 - SolvingProblem

  if (env.disorder_flg == false)
    ierr = VecDestroy(&real_vec); CHKERRQ(ierr);

#ifdef TIME_CODE
  env.tm.PrintTimeInfoFull();
#endif

  return ierr;
}
