#include "hamiltonian.hpp"
#include <slepcsys.h>

static char help[] = "This is a test to test the Hamiltonian<chain1D> class\n\n";

int main(int argc, char ** argv) {

  PetscReal J1 = 0.0;
  PetscReal Delta1 = 0.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;

  Environment env{argc,argv,6,help};
    
  Basis b{env,-1};
  chain1D lat{env};
  Hamiltonian<chain1D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;

  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  PetscInt * d_nnz, * o_nnz;
  PetscInt size = h.hamilt_size();
  PetscInt local_size = PETSC_DECIDE;
  ierr = PetscSplitOwnership(MPI_COMM_WORLD, &local_size, &size); CHKERRQ(ierr);
  PetscCalloc1(local_size, &d_nnz);
  PetscCalloc1(local_size, &o_nnz);

  for (PetscInt i = 0; i < local_size; ++i)
    d_nnz[i] = 1;
  ierr = MatMPIAIJSetPreallocation(h.hamilt, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);
  
  ierr = h.build_disorder(env,hi); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(MPI_COMM_WORLD, "build_disorder6m1.dat", &viewer);
  ierr = MatView(h.hamilt, viewer); CHKERRQ(ierr);

  ierr = PetscFree(d_nnz);
  ierr = PetscFree(o_nnz);
  ierr = PetscFree(hi);
  ierr = PetscViewerDestroy(&viewer);
  
  return ierr;
}
