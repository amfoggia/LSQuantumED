#include <parse_args.hpp>
#include <print_args.hpp>
#include <hamiltonian.hpp>
#include <random_disorder.hpp>
#include <solver.hpp>

static char help[] = "Diagonalization of XXZ-Heisenberg Hamiltonian with next-nearest neighbours interactions in a 1D chain lattice\n.";

int main(int argc, char * argv[]) {

  PetscErrorCode ierr;
  EPS solver;
  PetscInt nconv;
  double eigen, error;
  Vec state;

  // Initialize the environment
  Environment env{argc,argv,help};
  
  // Read command line arguments
  parse_args();
  Delta1 = 0.0;
  J1 = 0.0;
  J2 = 1.0;
  min_dis = -0.3;
  max_dis = 0.4;
  reprod = PETSC_FALSE;
  print_args();
  
  // Create Basis
  Basis b{env,0};

  // Create Lattice
  chain1D lat{env};

  // Create Hamiltonian
  Hamiltonian<chain1D> h{env,b,lat,J1,Delta1,J2,Delta2};
  ierr = h.build_diag(env); CHKERRQ(ierr);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);

  // Create random disorder strengths
  Tools::RandomDisorder rd{env.nspins,reprod,min_dis,max_dis,rep};

  ierr = MatCreateVecs(h.hamilt,&state,NULL); CHKERRQ(ierr);
  // Solve for every disorder realization
  for (int i = 0; i < rep; ++i) {
    rd.hi_init(i);
    ierr = h.build_disorder(env,rd.hi); CHKERRQ(ierr);
    ierr = Solver::SolverInit(solver,h.hamilt,10); CHKERRQ(ierr);
    ierr = Solver::solve(env,solver,nconv); CHKERRQ(ierr);
    ierr = Solver::solution(solver,0,eigen,state,error); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%d      %10.5f\n",i,eigen); CHKERRQ(ierr);
    ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
  }

  ierr = VecDestroy(&state); CHKERRQ(ierr);
  return ierr;
}

/*
Run with: -n <num_spins> -ltype chain1D -disorder -rep <num_repetitions> -d2 <d2_value>
*/
