#include "hamiltonian.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Hamiltonian class\n\n";
int _argc;
char ** _argv;

Environment env{_argc,_argv,6,help};

class HamiltonianTestEnv : public ::testing::Environment {
protected:
  virtual void TearDown() {
    ::testing::internal::CaptureStdout();
    EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};

TEST(PreallocInfo, chain1D) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 1.0;
    
    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 20;
    PetscInt local_size = PETSC_DECIDE;
    PetscInt rest = size%mpi_size;
    PetscInt gstart, gend;
    PetscSplitOwnership(MPI_COMM_WORLD, &local_size, &size);
    gstart = local_size*mpi_rank + (rest <= mpi_rank)*rest;
    gend = local_size*(mpi_rank+1) + (rest <= mpi_rank)*rest;

    PetscInt * d_nnz, * o_nnz, * coup_elems;
    PetscCalloc1(local_size, &d_nnz);
    PetscCalloc1(local_size, &o_nnz);
    PetscCalloc1(8, &coup_elems);
  
    // We are running with 3 processes (in meson.build)
    if (mpi_rank == 0) {

      coup_elems[0] = 1;
      coup_elems[1] = 12;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 2;
      coup_elems[1] = 0;
      coup_elems[2] = 4;
      coup_elems[3] = 14;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);

      coup_elems[0] = 3;
      coup_elems[1] = 1;
      coup_elems[2] = 5;
      coup_elems[3] = 15;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[2]);
      EXPECT_EQ(1,o_nnz[2]);

      coup_elems[0] = 2;
      coup_elems[1] = 6;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(0,o_nnz[3]);
    
      coup_elems[0] = 5;
      coup_elems[1] = 1;
      coup_elems[2] = 10;
      coup_elems[3] = 17;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(2,o_nnz[4]);

      coup_elems[0] = 6;
      coup_elems[1] = 4;
      coup_elems[2] = 7;
      coup_elems[3] = 2;
      coup_elems[4] = 11;
      coup_elems[5] = 18;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,6,-6,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[5]);
      EXPECT_EQ(3,o_nnz[5]);

      coup_elems[0] = 5;
      coup_elems[1] = 8;
      coup_elems[2] = 3;
      coup_elems[3] = 13;
      HamiltHelper::prealloc_info(Delta,6,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(2,o_nnz[4]);
    }

    if (mpi_rank == 1) {

      coup_elems[0] = 8;
      coup_elems[1] = 5;
      coup_elems[2] = 13;
      coup_elems[3] = 19;
      HamiltHelper::prealloc_info(Delta,7,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[0]);
      EXPECT_EQ(2,o_nnz[0]);

      coup_elems[0] = 7;
      coup_elems[1] = 9;
      coup_elems[2] = 6;
      coup_elems[3] = 14;
      HamiltHelper::prealloc_info(Delta,8,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 8;
      coup_elems[1] = 15;
      HamiltHelper::prealloc_info(Delta,9,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[2]);
      EXPECT_EQ(1,o_nnz[2]);

      coup_elems[0] = 11;
      coup_elems[1] = 4;
      HamiltHelper::prealloc_info(Delta,10,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[3]);
      EXPECT_EQ(1,o_nnz[3]);
    
      coup_elems[0] = 12;
      coup_elems[1] = 10;
      coup_elems[2] = 13;
      coup_elems[3] = 5;
      HamiltHelper::prealloc_info(Delta,11,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[4]);
      EXPECT_EQ(1,o_nnz[4]);

      coup_elems[0] = 11;
      coup_elems[1] = 14;
      coup_elems[2] = 6;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,12,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[5]);
      EXPECT_EQ(3,o_nnz[5]);

      coup_elems[0] = 14;
      coup_elems[1] = 11;
      coup_elems[2] = 16;
      coup_elems[3] = 7;
      HamiltHelper::prealloc_info(Delta,13,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[6]);
      EXPECT_EQ(2,o_nnz[6]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 13;
      coup_elems[1] = 15;
      coup_elems[2] = 12;
      coup_elems[3] = 17;
      coup_elems[4] = 8;
      coup_elems[5] = 1;
      HamiltHelper::prealloc_info(Delta,14,gstart,gend,coup_elems,6,-6,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[0]);
      EXPECT_EQ(4,o_nnz[0]);

      coup_elems[0] = 14;
      coup_elems[1] = 18;
      coup_elems[2] = 9;
      coup_elems[3] = 2;
      HamiltHelper::prealloc_info(Delta,15,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 17;
      coup_elems[1] = 13;
      HamiltHelper::prealloc_info(Delta,16,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[2]);
      EXPECT_EQ(1,o_nnz[2]);

      coup_elems[0] = 16;
      coup_elems[1] = 18;
      coup_elems[2] = 14;
      coup_elems[3] = 4;
      HamiltHelper::prealloc_info(Delta,17,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[3]);
      EXPECT_EQ(1,o_nnz[3]);
    
      coup_elems[0] = 17;
      coup_elems[1] = 19;
      coup_elems[2] = 15;
      coup_elems[3] = 5;
      HamiltHelper::prealloc_info(Delta,18,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[4]);
      EXPECT_EQ(1,o_nnz[4]);

      coup_elems[0] = 18;
      coup_elems[1] = 7;
      HamiltHelper::prealloc_info(Delta,19,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[5]);
      EXPECT_EQ(1,o_nnz[5]);
    }

    ierr = PetscFree(coup_elems); EXPECT_EQ(0,ierr);
    ierr = PetscFree(d_nnz); EXPECT_EQ(0,ierr);
    ierr = PetscFree(o_nnz); EXPECT_EQ(0,ierr);
  }
}

TEST(PreallocInfo, chain1D61) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 1.0;
    
    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 15;
    PetscInt local_size = PETSC_DECIDE;
    PetscInt rest = size%mpi_size;
    PetscInt gstart, gend;
    PetscSplitOwnership(MPI_COMM_WORLD, &local_size, &size);
    gstart = local_size*mpi_rank + (rest <= mpi_rank)*rest;
    gend = local_size*(mpi_rank+1) + (rest <= mpi_rank)*rest;

    PetscInt * d_nnz, * o_nnz, * coup_elems;
    PetscCalloc1(local_size, &d_nnz);
    PetscCalloc1(local_size, &o_nnz);
    PetscCalloc1(8, &coup_elems);
  
    // We are running with 3 processes (in meson.build)
    if (mpi_rank == 0) {

      coup_elems[0] = 1;
      coup_elems[1] = 8;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 2;
      coup_elems[1] = 0;
      coup_elems[2] = 5;
      coup_elems[3] = 11;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 3;
      coup_elems[1] = 1;
      coup_elems[2] = 6;
      coup_elems[3] = 13;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(2,o_nnz[2]);

      coup_elems[0] = 4;
      coup_elems[1] = 2;
      coup_elems[2] = 7;
      coup_elems[3] = 14;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);
    
      coup_elems[0] = 3;
      coup_elems[1] = 8;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[4]);
      EXPECT_EQ(1,o_nnz[4]);
    }

    if (mpi_rank == 1) {

      coup_elems[0] = 6;
      coup_elems[1] = 1;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 7;
      coup_elems[1] = 5;
      coup_elems[2] = 9;
      coup_elems[3] = 2;
      HamiltHelper::prealloc_info(Delta,6,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);

      coup_elems[0] = 8;
      coup_elems[1] = 6;
      coup_elems[2] = 10;
      coup_elems[3] = 3;
      HamiltHelper::prealloc_info(Delta,7,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(2,o_nnz[2]);

      coup_elems[0] = 7;
      coup_elems[1] = 11;
      coup_elems[2] = 4;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,8,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[3]);
      EXPECT_EQ(3,o_nnz[3]);
    
      coup_elems[0] = 10;
      coup_elems[1] = 6;
      HamiltHelper::prealloc_info(Delta,9,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[4]);
      EXPECT_EQ(1,o_nnz[4]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 11;
      coup_elems[1] = 9;
      coup_elems[2] = 12;
      coup_elems[3] = 7;
      HamiltHelper::prealloc_info(Delta,10,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[0]);
      EXPECT_EQ(2,o_nnz[0]);

      coup_elems[0] = 10;
      coup_elems[1] = 13;
      coup_elems[2] = 8;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,11,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 13;
      coup_elems[1] = 10;
      HamiltHelper::prealloc_info(Delta,12,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(0,o_nnz[2]);

      coup_elems[0] = 12;
      coup_elems[1] = 14;
      coup_elems[2] = 11;
      coup_elems[3] = 2;
      HamiltHelper::prealloc_info(Delta,13,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[3]);
      EXPECT_EQ(1,o_nnz[3]);
    
      coup_elems[0] = 13;
      coup_elems[1] = 3;
      HamiltHelper::prealloc_info(Delta,14,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[4]);
      EXPECT_EQ(1,o_nnz[4]);
    }

    ierr = PetscFree(coup_elems); EXPECT_EQ(0,ierr);
    ierr = PetscFree(d_nnz); EXPECT_EQ(0,ierr);
    ierr = PetscFree(o_nnz); EXPECT_EQ(0,ierr);
  }
}

TEST(PreallocInfo, chain1D6m1) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 1.0;
    
    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 15;
    PetscInt local_size = PETSC_DECIDE;
    PetscInt rest = size%mpi_size;
    PetscInt gstart, gend;
    PetscSplitOwnership(MPI_COMM_WORLD, &local_size, &size);
    gstart = local_size*mpi_rank + (rest <= mpi_rank)*rest;
    gend = local_size*(mpi_rank+1) + (rest <= mpi_rank)*rest;

    PetscInt * d_nnz, * o_nnz, * coup_elems;
    PetscCalloc1(local_size, &d_nnz);
    PetscCalloc1(local_size, &o_nnz);
    PetscCalloc1(8, &coup_elems);
  
    // We are running with 3 processes (in meson.build)
    if (mpi_rank == 0) {

      coup_elems[0] = 1;
      coup_elems[1] = 11;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 2;
      coup_elems[1] = 0;
      coup_elems[2] = 3;
      coup_elems[3] = 12;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);

      coup_elems[0] = 1;
      coup_elems[1] = 4;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(0,o_nnz[2]);

      coup_elems[0] = 4;
      coup_elems[1] = 1;
      coup_elems[2] = 6;
      coup_elems[3] = 13;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);
    
      coup_elems[0] = 3;
      coup_elems[1] = 5;
      coup_elems[2] = 2;
      coup_elems[3] = 7;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(2,o_nnz[4]);
    }

    if (mpi_rank == 1) {

      coup_elems[0] = 4;
      coup_elems[1] = 8;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 7;
      coup_elems[1] = 3;
      coup_elems[2] = 10;
      coup_elems[3] = 14;
      HamiltHelper::prealloc_info(Delta,6,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[1]);
      EXPECT_EQ(3,o_nnz[1]);

      coup_elems[0] = 6;
      coup_elems[1] = 8;
      coup_elems[2] = 4;
      coup_elems[3] = 11;
      HamiltHelper::prealloc_info(Delta,7,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(2,o_nnz[2]);

      coup_elems[0] = 7;
      coup_elems[1] = 9;
      coup_elems[2] = 5;
      coup_elems[3] = 12;
      HamiltHelper::prealloc_info(Delta,8,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[3]);
      EXPECT_EQ(1,o_nnz[3]);
    
      coup_elems[0] = 8;
      coup_elems[1] = 13;
      HamiltHelper::prealloc_info(Delta,9,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[4]);
      EXPECT_EQ(1,o_nnz[4]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 11;
      coup_elems[1] = 6;
      HamiltHelper::prealloc_info(Delta,10,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 10;
      coup_elems[1] = 12;
      coup_elems[2] = 7;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,11,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 11;
      coup_elems[1] = 13;
      coup_elems[2] = 8;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,12,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(2,o_nnz[2]);

      coup_elems[0] = 12;
      coup_elems[1] = 14;
      coup_elems[2] = 9;
      coup_elems[3] = 3;
      HamiltHelper::prealloc_info(Delta,13,gstart,gend,coup_elems,4,-2,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);
    
      coup_elems[0] = 13;
      coup_elems[1] = 6;
      HamiltHelper::prealloc_info(Delta,14,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[4]);
      EXPECT_EQ(1,o_nnz[4]);
    }

    ierr = PetscFree(coup_elems); EXPECT_EQ(0,ierr);
    ierr = PetscFree(d_nnz); EXPECT_EQ(0,ierr);
    ierr = PetscFree(o_nnz); EXPECT_EQ(0,ierr);
  }
}

TEST(PreallocInfo, chain1D62) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 1.0;
    
    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 6;
    PetscInt local_size = PETSC_DECIDE;
    PetscInt rest = size%mpi_size;
    PetscInt gstart, gend;
    PetscSplitOwnership(MPI_COMM_WORLD, &local_size, &size);
    gstart = local_size*mpi_rank + (rest <= mpi_rank)*rest;
    gend = local_size*(mpi_rank+1) + (rest <= mpi_rank)*rest;

    PetscInt * d_nnz, * o_nnz, * coup_elems;
    PetscCalloc1(local_size, &d_nnz);
    PetscCalloc1(local_size, &o_nnz);
    PetscCalloc1(8, &coup_elems);
  
    // We are running with 3 processes (in meson.build)
    if (mpi_rank == 0) {

      coup_elems[0] = 1;
      coup_elems[1] = 5;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 2;
      coup_elems[1] = 0;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);
    }

    if (mpi_rank == 1) {

      coup_elems[0] = 3;
      coup_elems[1] = 1;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 4;
      coup_elems[1] = 2;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 5;
      coup_elems[1] = 3;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 4;
      coup_elems[1] = 0;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);
    }

    ierr = PetscFree(coup_elems); EXPECT_EQ(0,ierr);
    ierr = PetscFree(d_nnz); EXPECT_EQ(0,ierr);
    ierr = PetscFree(o_nnz); EXPECT_EQ(0,ierr);
  }
}

TEST(PreallocInfo, chain1D6m2) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 1.0;
    
    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 6;
    PetscInt local_size = PETSC_DECIDE;
    PetscInt rest = size%mpi_size;
    PetscInt gstart, gend;
    PetscSplitOwnership(MPI_COMM_WORLD, &local_size, &size);
    gstart = local_size*mpi_rank + (rest <= mpi_rank)*rest;
    gend = local_size*(mpi_rank+1) + (rest <= mpi_rank)*rest;

    PetscInt * d_nnz, * o_nnz, * coup_elems;
    PetscCalloc1(local_size, &d_nnz);
    PetscCalloc1(local_size, &o_nnz);
    PetscCalloc1(8, &coup_elems);
  
    // We are running with 3 processes (in meson.build)
    if (mpi_rank == 0) {

      coup_elems[0] = 1;
      coup_elems[1] = 5;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 0;
      coup_elems[1] = 2;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);
    }

    if (mpi_rank == 1) {

      coup_elems[0] = 1;
      coup_elems[1] = 3;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 2;
      coup_elems[1] = 4;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 3;
      coup_elems[1] = 5;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(1,o_nnz[0]);

      coup_elems[0] = 4;
      coup_elems[1] = 0;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,2,2,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[1]);
      EXPECT_EQ(1,o_nnz[1]);
    }

    ierr = PetscFree(coup_elems); EXPECT_EQ(0,ierr);
    ierr = PetscFree(d_nnz); EXPECT_EQ(0,ierr);
    ierr = PetscFree(o_nnz); EXPECT_EQ(0,ierr);
  }
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

