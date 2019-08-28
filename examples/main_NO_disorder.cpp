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

  PetscErrorCode ierr = 0;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;
  EPS solver;
  Vec real_vec;
    
  Environment env{argc,argv,help};
  
  parse_args();
  print_args();
  
  Basis b{env,0};
  chain1D lattice{env};
  Hamiltonian<chain1D> h{env,b,lattice,J1,Delta1,J2,Delta2};
  ierr = h.build_diag(env); CHKERRQ(ierr);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);
  ierr = Solver::SolverInit(solver,h.hamilt,1); CHKERRQ(ierr);
  ierr = MatCreateVecs(h.hamilt,&real_vec,NULL); CHKERRQ(ierr);

  ierr = Solver::solve(env,solver,nconv); CHKERRQ(ierr);

  if (nconv > 0) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
		       "           k          ||Ax-kx||/||kx||\n"
		       "   ----------------- ------------------\n"); CHKERRQ(ierr);
    
    for (PetscInt i = 0; i < 1; ++i) {
      ierr = Solver::solution(solver,i,real_eig,real_vec,error); CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",PetscRealPart(real_eig),(double)error); CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);
  }

  ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
  ierr = VecDestroy(&real_vec); CHKERRQ(ierr);

  return ierr;
}
