#include "hamiltonian.hpp"
#include <slepcsys.h>

static char help[] = "This is a test to test the Hamiltonian<chain1D> class\n\n";

int main(int argc, char ** argv) {

  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  
  Environment env{argc,argv,6,help};
    
  Basis b{env,3};
  chain1D lat{env};
  Hamiltonian<chain1D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;

  ierr = h.build_off_diag(env); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(MPI_COMM_WORLD, "build_off_diag63.dat", &viewer);
  ierr = MatView(h.hamilt, viewer); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  return ierr;
}
