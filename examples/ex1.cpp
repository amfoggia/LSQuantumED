#include <parse_args.hpp>
#include <print_args.hpp>
#include <hamiltonian.hpp>
#include <solver.hpp>

static char help[] = "Diagonalization of XXZ-Heisenberg Hamiltonian with nearest neighbours interactions in a 2D square lattice\n.";

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
  Delta1 = 0.5;
  J2 = 0.0;
  Delta2 = 0.0;
  print_args();
  
  // Create Basis
  Basis b{env,0};

  // Create Lattice
  square2D lat{env,lx,ly};

  // Create Hamiltonian
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  ierr = h.build_diag(env); CHKERRQ(ierr);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);
  // ierr = MatView(h.hamilt, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

  // Initiate solver and solve
  ierr = Solver::SolverInit(solver,h.hamilt,10); CHKERRQ(ierr);
  ierr = Solver::solve(env,solver,nconv); CHKERRQ(ierr);

  // Print eigenvalues
  ierr = MatCreateVecs(h.hamilt,&state,NULL); CHKERRQ(ierr);
  for (int i = 0; i < 10; ++i) {
    ierr = Solver::solution(solver,i,eigen,state,error); CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%d      %10.5f\n",i,eigen); CHKERRQ(ierr);
  }

  ierr = VecDestroy(&state); CHKERRQ(ierr);
  ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
  return ierr;
}

/*
Run with: -n <num_spins> -ltype square2D -lx <x-dim> -ly <y-dim>
*/
