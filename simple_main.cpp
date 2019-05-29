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

  Tools::RandomDisorder rd{spins,reprod,min_dis,max_dis,rep};
  
  Basis b{env,0};

  chain1D lattice{env};
  Hamiltonian<chain1D> h{env,b,lattice,J1,Delta1,J2,Delta2};
  ierr = h.build_diag(env); CHKERRQ(ierr);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);

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
  dsf_data.nev = 1;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
  dsf_data.disorder = PETSC_TRUE;
  dsf_data.path = "./";

  for (int repet = 0; repet < rep; ++repet) {
    rd.hi_init(repet);
    ierr = h.build_disorder(env,rd.hi); CHKERRQ(ierr);
    dsf_data.hi = rd.hi;
    dsf_data.h0 = &h;
    dsf_data.h1 = &h;
    ierr = Solver::SolverInit(solver,h.hamilt,1); CHKERRQ(ierr);
    ierr = Solver::solve_lanczos(env,solver,nconv); CHKERRQ(ierr);
    ierr = MatCreateVecs(h.hamilt,&real_vec,NULL); CHKERRQ(ierr);
    if (nconv > 0) {      
      for (PetscInt i = 0; i < 1; ++i) {
	ierr = Solver::solution(solver,i,real_eig,real_vec,error); CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"  %d     %12f       %12g\n",repet,PetscRealPart(real_eig),(double)error);
	CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);
    } // -- if nconv

    dynstructfactor = Phys::DynStructFactor<chain1D,Sqz>(env, dsf_data, real_eig, real_vec, repet);

    ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
    ierr = VecDestroy(&real_vec); CHKERRQ(ierr);
  } // -- for repet

  return ierr;
}
