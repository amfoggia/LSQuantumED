#include "hamiltonian.hpp"
#include <slepcsys.h>

static char help[] = "This is a test to test the Hamiltonian<square2D> class\n\n";

int main(int argc, char ** argv) {

  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  
  Environment env{argc,argv,12,help};
    
  Basis b{env,5};
  square2D lat{env,3,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  
  ierr = h.build_off_diag(env); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(MPI_COMM_WORLD, "build_off_diag125.dat", &viewer);
  ierr = MatView(h.hamilt, viewer); CHKERRQ(ierr);

  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  return ierr;
}
