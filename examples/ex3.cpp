#include <parse_args.hpp>
#include <print_args.hpp>
#include <hamiltonian.hpp>
#include <solver.hpp>
#include <phys.hpp>

static char help[] = "Diagonalization of XXZ-Heisenberg Hamiltonian with nearest neighbours interactions in a 2D square lattice\n.Computation of physical quantities\n.";

int main(int argc, char * argv[]) {

  PetscErrorCode ierr;
  EPS solver;
  PetscInt nconv;
  double eigen, error;
  Vec state;
  double magAF, magSTR;

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

  // Create sublattices
  AF<square2D> sub_af{lat,2};
  Striped<square2D> sub_str{lat,4};

  // Create Hamiltonian
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  ierr = h.build_diag(env); CHKERRQ(ierr);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);

  // Initiate solver and solve
  ierr = Solver::SolverInit(solver,h.hamilt,1); CHKERRQ(ierr);
  ierr = Solver::solve(env,solver,nconv); CHKERRQ(ierr);

  // Print groundstate
  ierr = MatCreateVecs(h.hamilt,&state,NULL); CHKERRQ(ierr);
  ierr = Solver::solution(solver,0,eigen,state,error); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%d      %10.5f\n",0,eigen); CHKERRQ(ierr);

  // Compute order parameters
  magAF = Phys::SquaredSubLattMagAF<square2D>(env,b,sub_af,state);
  magSTR = Phys::SquaredSubLattMagSTR<square2D>(env,b,sub_str,state);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"magAF: %.10f\n",magAF); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"magSTR: %.10f\n",magSTR); CHKERRQ(ierr);

  // Create data object for DSF
  Phys::DSF_data<square2D> dsf_data{};
  dsf_data.b0 = &b;
  dsf_data.b1 = &b;
  dsf_data.lat = &lat;
  dsf_data.h0 = &h;
  dsf_data.h1 = &h;
  dsf_data.nev = 10;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
  dsf_data.disorder = PETSC_FALSE;
  dsf_data.hi = NULL;
  dsf_data.path = "./";

  // Compute DSF
  ierr = Phys::DynStructFactor<square2D,Sqz>(env,dsf_data,eigen,state,0);
  
  ierr = VecDestroy(&state); CHKERRQ(ierr);
  ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
  return ierr;
}

/*
Run with: -n <num_spins> -ltype square2D -lx <x-dim> -ly <y-dim>
*/
