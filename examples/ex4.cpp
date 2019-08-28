#include <parse_args.hpp>
#include <print_args.hpp>
#include "solver.hpp"
#include "phys.hpp"

static char help[] = "Computation of dynamical structure factor with S+ operator.\n";

int main (int argc, char * argv[]) {

  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;
  
  Environment env{argc,argv,help};

  parse_args();
  print_args();
  
  Basis b{env,0};
  square2D lat{env,lx,ly};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};  

  ierr = h.build_diag(env); CHKERRQ(ierr);
  ierr = h.build_off_diag(env); CHKERRQ(ierr);
  
  ierr = Solver::SolverInit(solver,h.hamilt,1); CHKERRQ(ierr);
  ierr = Solver::solve(env,solver,nconv); CHKERRQ(ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); CHKERRQ(ierr);

  ierr = Solver::solution(solver, 0, real_eig, xr, error); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "real_eig: %f\n", real_eig);

  Basis b1{env,1};
  Hamiltonian<square2D> h1{env,b1,lat,J1,Delta1,J2,Delta2};
  ierr = h1.build_diag(env); CHKERRQ(ierr);
  ierr = h1.build_off_diag(env); CHKERRQ(ierr);

  Phys::DSF_data<square2D> dsf_data{};
  dsf_data.b0 = &b;
  dsf_data.b1 = &b1;
  dsf_data.lat = &lat;
  dsf_data.nev = 10;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
  dsf_data.h0 = &h;
  dsf_data.h1 = &h1;
  dsf_data.disorder = PETSC_FALSE;
  dsf_data.hi = NULL;
  dsf_data.path = "./";
  
  ierr = Phys::DynStructFactor<square2D,Sqp>(env, dsf_data, real_eig, xr, 0); CHKERRQ(ierr);
  
  ierr = Solver::SolverClean(solver); CHKERRQ(ierr);
  ierr = VecDestroy(&xr); CHKERRQ(ierr);

  return ierr;
  
}
