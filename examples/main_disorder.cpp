#include "sublattice.hpp"
#include "time_monitor.hpp"
#include "random_disorder.hpp"
#include "parse_args.hpp"
#include "print_args.hpp"
#include "phys.hpp"
#include <cstdlib>
#include <petsctime.h>

static char help[] = "Hamiltonian with constant total magnetization (U(1) symmetry)\n";

int main(int argc, char * argv[]) {

  EPS solver;

  PetscErrorCode ierr = 0;
  PetscLogStage stage1, stage2, stage3, stage4a, stage4b;
  PetscInt nconv, spins;
  PetscScalar real_eig;
  Vec real_vec;
  PetscReal error;
  PetscScalar exp_val, magAF = 0.0;
  
  Environment env{argc,argv,spins,help};
  
  parse_args();
  print_args();
  
  Tools::RandomDisorder rd{spins,reprod,min_dis,max_dis,rep};

  ierr = PetscPrintf(MPI_COMM_WORLD, "1- Creating basis\n");
  
  PetscLogStageRegister("Create Basis", &stage1);
  PetscLogStagePush(stage1);
  Basis b{env,0};
  PetscLogStagePop();
  
  ierr = PetscPrintf(MPI_COMM_WORLD, "2- Creating lattice\n"); CHKERRQ(ierr);

  PetscLogStageRegister("Create Latice", &stage2);
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
    
  //  ierr = PetscPrintf(MPI_COMM_WORLD, "4c- Building disorder part\n"); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
		     " iter          k          ||Ax-kx||/||kx||\n"
		     " ----    -------------   ------------------\n"); CHKERRQ(ierr);

  PetscScalar dynstructfactor;
  Phys::DSF_data<chain1D> dsf_data{};
  dsf_data.b0 = &b;
  dsf_data.b1 = &b;
  dsf_data.lat = &lattice;
  dsf_data.nev = 10;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
  dsf_data.disorder = PETSC_TRUE;
  dsf_data.path = "./";
  
  for (int repet = 0; repet < rep; ++repet) {
    rd.hi_init(repet);
    ierr = h.build_disorder(env,rd.hi); CHKERRQ(ierr);
    dsf_data.h0 = &h;
    dsf_data.h1 = &h;
    dsf_data.hi = rd.hi;
    ierr = Solver::SolverInit(solver,h.hamilt,1); CHKERRQ(ierr);
    ierr = Solver::solve_lanczos(env,solver,nconv); CHKERRQ(ierr);
    ierr = MatCreateVecs(h.hamilt,&real_vec,NULL); CHKERRQ(ierr); // --> In here, only when there's disorder
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

    dynstructfactor = Phys::DynStructFactor<chain1D,Sqz>(env, dsf_data, real_eig, real_vec, repet);

    ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
    ierr = VecDestroy(&real_vec); CHKERRQ(ierr); // --> In here, only when there's disorder
  } // -- for repet

  ierr = PetscPrintf(MPI_COMM_WORLD, "9- Cleaning\n"); CHKERRQ(ierr);

#ifdef TIME_CODE
  env.tm.PrintTimeInfoFull();
#endif

  return ierr;
}
