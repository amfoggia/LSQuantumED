#include "hamiltonian.hpp"
#include <slepcsys.h>

static char help[] = "This is a test to test the Hamiltonian class\n\n";

int main(int argc, char ** argv) {

  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 1.0;
  PetscReal Delta2 = 1.0;

  Environment env{argc,argv,12,help};
    
  Basis b{env,5};
  square2D lat{env,3,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  env.disorder_flg = PETSC_TRUE;

  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  ierr = h.build_off_diag(env); CHKERRQ(ierr);
  ierr = h.build_disorder(env,hi); CHKERRQ(ierr);
  
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(MPI_COMM_WORLD, "build_disorder125.dat", &viewer);
  ierr = MatView(h.hamilt, viewer); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
  ierr = PetscFree(hi); CHKERRQ(ierr);

  return ierr;
}
