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

Environment env{_argc,_argv,8,help};

class HamiltonianTestEnv : public ::testing::Environment {
protected:
  virtual void TearDown() {
    ::testing::internal::CaptureStdout();
    EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};

TEST(PreallocInfo, square2D82) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 1.0;
    
    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 28;
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
  
    // We are running with 5 processes (in meson.build)
    if (mpi_rank == 0) {

      coup_elems[0] = 4;
      coup_elems[1] = 9;
      coup_elems[2] = 1;
      coup_elems[3] = 4;
      coup_elems[4] = 8;
      coup_elems[5] = 9;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[0]);
      EXPECT_EQ(3,o_nnz[0]);

      coup_elems[0] = 5;
      coup_elems[1] = 14;
      coup_elems[2] = 2;
      coup_elems[3] = 0;
      coup_elems[4] = 5;
      coup_elems[5] = 7;
      coup_elems[6] = 13;
      coup_elems[7] = 14;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[1]);
      EXPECT_EQ(4,o_nnz[1]);

      coup_elems[0] = 6;
      coup_elems[1] = 18;
      coup_elems[2] = 1;
      coup_elems[3] = 6;
      coup_elems[4] = 8;
      coup_elems[5] = 18;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[2]);
      EXPECT_EQ(5,o_nnz[2]);

      coup_elems[0] = 4;
      coup_elems[1] = 6;
      coup_elems[2] = 9;
      coup_elems[3] = 18;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[3]);
      EXPECT_EQ(3,o_nnz[3]);
    
      coup_elems[0] = 5;
      coup_elems[1] = 3;
      coup_elems[2] = 0;
      coup_elems[3] = 22;
      coup_elems[4] = 10;
      coup_elems[5] = 0;
      coup_elems[6] = 19;
      coup_elems[7] = 22;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[4]);
      EXPECT_EQ(4,o_nnz[4]);

      coup_elems[0] = 6;
      coup_elems[1] = 4;
      coup_elems[2] = 1;
      coup_elems[3] = 23;
      coup_elems[4] = 1;
      coup_elems[5] = 11;
      coup_elems[6] = 20;
      coup_elems[7] = 23;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[5]);
      EXPECT_EQ(5,o_nnz[5]);
    }

    if (mpi_rank == 1) {

      coup_elems[0] = 5;
      coup_elems[1] = 2;
      coup_elems[2] = 3;
      coup_elems[3] = 24;
      coup_elems[4] = 2;
      coup_elems[5] = 12;
      coup_elems[6] = 21;
      coup_elems[7] = 24;
      HamiltHelper::prealloc_info(Delta,6,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(8,o_nnz[0]);

      coup_elems[0] = 11;
      coup_elems[1] = 15;
      coup_elems[2] = 8;
      coup_elems[3] = 11;
      coup_elems[4] = 1;
      coup_elems[5] = 15;
      HamiltHelper::prealloc_info(Delta,7,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[1]);
      EXPECT_EQ(3,o_nnz[1]);

      coup_elems[0] = 12;
      coup_elems[1] = 19;
      coup_elems[2] = 7;
      coup_elems[3] = 12;
      coup_elems[4] = 13;
      coup_elems[5] = 2;
      coup_elems[6] = 19;
      coup_elems[7] = 0;
      HamiltHelper::prealloc_info(Delta,8,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[2]);
      EXPECT_EQ(7,o_nnz[2]);

      coup_elems[0] = 10;
      coup_elems[1] = 22;
      coup_elems[2] = 12;
      coup_elems[3] = 0;
      coup_elems[4] = 14;
      coup_elems[5] = 3;
      coup_elems[6] = 22;
      coup_elems[7] = 0;
      HamiltHelper::prealloc_info(Delta,9,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[3]);
      EXPECT_EQ(7,o_nnz[3]);
    
      coup_elems[0] = 11;
      coup_elems[1] = 9;
      coup_elems[2] = 15;
      coup_elems[3] = 4;
      HamiltHelper::prealloc_info(Delta,10,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(2,o_nnz[4]);

      coup_elems[0] = 12;
      coup_elems[1] = 10;
      coup_elems[2] = 7;
      coup_elems[3] = 25;
      coup_elems[4] = 16;
      coup_elems[5] = 7;
      coup_elems[6] = 5;
      coup_elems[7] = 25;
      HamiltHelper::prealloc_info(Delta,11,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[5]);
      EXPECT_EQ(5,o_nnz[5]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 11;
      coup_elems[1] = 8;
      coup_elems[2] = 26;
      coup_elems[3] = 9;
      coup_elems[4] = 8;
      coup_elems[5] = 17;
      coup_elems[6] = 6;
      coup_elems[7] = 26;
      HamiltHelper::prealloc_info(Delta,12,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(1,d_nnz[0]);
      EXPECT_EQ(7,o_nnz[0]);

      coup_elems[0] = 17;
      coup_elems[1] = 20;
      coup_elems[2] = 17;
      coup_elems[3] = 8;
      coup_elems[4] = 20;
      coup_elems[5] = 1;
      HamiltHelper::prealloc_info(Delta,13,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(4,o_nnz[1]);

      coup_elems[0] = 23;
      coup_elems[1] = 15;
      coup_elems[2] = 17;
      coup_elems[3] = 1;
      coup_elems[4] = 18;
      coup_elems[5] = 9;
      coup_elems[6] = 23;
      coup_elems[7] = 1;
      HamiltHelper::prealloc_info(Delta,14,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(6,o_nnz[2]);

      coup_elems[0] = 16;
      coup_elems[1] = 25;
      coup_elems[2] = 14;
      coup_elems[3] = 7;
      coup_elems[4] = 19;
      coup_elems[5] = 10;
      coup_elems[6] = 25;
      coup_elems[7] = 7;
      HamiltHelper::prealloc_info(Delta,15,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(6,o_nnz[3]);
    
      coup_elems[0] = 17;
      coup_elems[1] = 15;
      coup_elems[2] = 20;
      coup_elems[3] = 11;
      HamiltHelper::prealloc_info(Delta,16,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(2,o_nnz[4]);

      coup_elems[0] = 16;
      coup_elems[1] = 13;
      coup_elems[2] = 27;
      coup_elems[3] = 14;
      coup_elems[4] = 21;
      coup_elems[5] = 13;
      coup_elems[6] = 12;
      coup_elems[7] = 27;
      HamiltHelper::prealloc_info(Delta,17,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(5,d_nnz[5]);
      EXPECT_EQ(3,o_nnz[5]);
    }

    if (mpi_rank == 3) {

      coup_elems[0] = 24;
      coup_elems[1] = 19;
      coup_elems[2] = 21;
      coup_elems[3] = 2;
      coup_elems[4] = 14;
      coup_elems[5] = 24;
      coup_elems[6] = 3;
      coup_elems[7] = 2;
      HamiltHelper::prealloc_info(Delta,18,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[0]);
      EXPECT_EQ(6,o_nnz[0]);
      
      coup_elems[0] = 26;
      coup_elems[1] = 20;
      coup_elems[2] = 18;
      coup_elems[3] = 8;
      coup_elems[4] = 15;
      coup_elems[5] = 26;
      coup_elems[6] = 8;
      coup_elems[7] = 4;
      HamiltHelper::prealloc_info(Delta,19,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(6,o_nnz[1]);

      coup_elems[0] = 21;
      coup_elems[1] = 27;
      coup_elems[2] = 19;
      coup_elems[3] = 13;
      coup_elems[4] = 16;
      coup_elems[5] = 27;
      coup_elems[6] = 13;
      coup_elems[7] = 5;
      HamiltHelper::prealloc_info(Delta,20,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(6,o_nnz[2]);

      coup_elems[0] = 20;
      coup_elems[1] = 18;
      coup_elems[2] = 17;
      coup_elems[3] = 6;
      HamiltHelper::prealloc_info(Delta,21,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);

      coup_elems[0] = 23;
      coup_elems[1] = 9;
      coup_elems[2] = 26;
      coup_elems[3] = 4;
      coup_elems[4] = 9;
      coup_elems[5] = 4;
      HamiltHelper::prealloc_info(Delta,22,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[4]);
      EXPECT_EQ(6,o_nnz[4]);
    }

    if (mpi_rank == 4) {

      coup_elems[0] = 24;
      coup_elems[1] = 22;
      coup_elems[2] = 14;
      coup_elems[3] = 25;
      coup_elems[4] = 27;
      coup_elems[5] = 5;
      coup_elems[6] = 14;
      coup_elems[7] = 5;
      HamiltHelper::prealloc_info(Delta,23,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[0]);
      EXPECT_EQ(5,o_nnz[0]);

      coup_elems[0] = 23;
      coup_elems[1] = 18;
      coup_elems[2] = 26;
      coup_elems[3] = 6;
      coup_elems[4] = 18;
      coup_elems[5] = 6;
      HamiltHelper::prealloc_info(Delta,24,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(4,o_nnz[1]);

      coup_elems[0] = 26;
      coup_elems[1] = 15;
      coup_elems[2] = 23;
      coup_elems[3] = 11;
      coup_elems[4] = 15;
      coup_elems[5] = 11;
      HamiltHelper::prealloc_info(Delta,25,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(4,o_nnz[2]);

      coup_elems[0] = 25;
      coup_elems[1] = 19;
      coup_elems[2] = 27;
      coup_elems[3] = 24;
      coup_elems[4] = 12;
      coup_elems[5] = 22;
      coup_elems[6] = 19;
      coup_elems[7] = 12;
      HamiltHelper::prealloc_info(Delta,26,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(3,d_nnz[3]);
      EXPECT_EQ(5,o_nnz[3]);

      coup_elems[0] = 20;
      coup_elems[1] = 26;
      coup_elems[2] = 17;
      coup_elems[3] = 23;
      coup_elems[4] = 20;
      coup_elems[5] = 17;
      HamiltHelper::prealloc_info(Delta,27,gstart,gend,coup_elems,6,4,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(4,o_nnz[4]);
    }

    ierr = PetscFree(coup_elems); EXPECT_EQ(0,ierr);
    ierr = PetscFree(d_nnz); EXPECT_EQ(0,ierr);
    ierr = PetscFree(o_nnz); EXPECT_EQ(0,ierr);
  }
}

TEST(PreallocInfo, square2D82_J2) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 0.0;
    
    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 28;
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
  
    // We are running with 5 processes (in meson.build)
    if (mpi_rank == 0) {

      coup_elems[0] = 12;
      coup_elems[1] = 5;
      coup_elems[2] = 10;
      coup_elems[3] = 3;
      coup_elems[4] = 3;
      coup_elems[5] = 5;
      coup_elems[6] = 12;
      coup_elems[7] = 10;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[0]);
      EXPECT_EQ(4,o_nnz[0]);

      coup_elems[0] = 6;
      coup_elems[1] = 17;
      coup_elems[2] = 15;
      coup_elems[3] = 4;
      coup_elems[4] = 4;
      coup_elems[5] = 6;
      coup_elems[6] = 17;
      coup_elems[7] = 15;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(6,o_nnz[1]);

      coup_elems[0] = 21;
      coup_elems[1] = 5;
      coup_elems[2] = 19;
      coup_elems[3] = 3;
      coup_elems[4] = 5;
      coup_elems[5] = 3;
      coup_elems[6] = 21;
      coup_elems[7] = 19;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[2]);
      EXPECT_EQ(4,o_nnz[2]);

      coup_elems[0] = 24;
      coup_elems[1] = 22;
      coup_elems[2] = 2;
      coup_elems[3] = 0;
      coup_elems[4] = 2;
      coup_elems[5] = 0;
      coup_elems[6] = 24;
      coup_elems[7] = 22;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[3]);
      EXPECT_EQ(4,o_nnz[3]);
    
      coup_elems[0] = 26;
      coup_elems[1] = 1;
      coup_elems[2] = 1;
      coup_elems[3] = 26;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(2,o_nnz[4]);

      coup_elems[0] = 27;
      coup_elems[1] = 0;
      coup_elems[2] = 2;
      coup_elems[3] = 25;
      coup_elems[4] = 2;
      coup_elems[5] = 0;
      coup_elems[6] = 27;
      coup_elems[7] = 25;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[5]);
      EXPECT_EQ(4,o_nnz[5]);
    }

    if (mpi_rank == 1) {

      coup_elems[0] = 1;
      coup_elems[1] = 26;
      coup_elems[2] = 1;
      coup_elems[3] = 26;
      HamiltHelper::prealloc_info(Delta,6,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(4,o_nnz[0]);

      coup_elems[0] = 12;
      coup_elems[1] = 16;
      coup_elems[2] = 10;
      coup_elems[3] = 14;
      coup_elems[4] = 10;
      coup_elems[5] = 12;
      coup_elems[6] = 14;
      coup_elems[7] = 16;
      HamiltHelper::prealloc_info(Delta,7,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(6,o_nnz[1]);

      coup_elems[0] = 20;
      coup_elems[1] = 11;
      coup_elems[2] = 9;
      coup_elems[3] = 18;
      coup_elems[4] = 11;
      coup_elems[5] = 9;
      coup_elems[6] = 18;
      coup_elems[7] = 20;
      HamiltHelper::prealloc_info(Delta,8,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[2]);
      EXPECT_EQ(4,o_nnz[2]);

      coup_elems[0] = 23;
      coup_elems[1] = 8;
      coup_elems[2] = 8;
      coup_elems[3] = 23;
      HamiltHelper::prealloc_info(Delta,9,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);
    
      coup_elems[0] = 25;
      coup_elems[1] = 0;
      coup_elems[2] = 7;
      coup_elems[3] = 22;
      coup_elems[4] = 7;
      coup_elems[5] = 22;
      coup_elems[6] = 25;
      coup_elems[7] = 0;
      HamiltHelper::prealloc_info(Delta,10,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(6,o_nnz[4]);

      coup_elems[0] = 8;
      coup_elems[1] = 23;
      coup_elems[2] = 8;
      coup_elems[3] = 23;
      HamiltHelper::prealloc_info(Delta,11,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[5]);
      EXPECT_EQ(2,o_nnz[5]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 7;
      coup_elems[1] = 0;
      coup_elems[2] = 27;
      coup_elems[3] = 24;
      coup_elems[4] = 7;
      coup_elems[5] = 24;
      coup_elems[6] = 27;
      coup_elems[7] = 0;
      HamiltHelper::prealloc_info(Delta,12,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(8,o_nnz[0]);

      coup_elems[0] = 21;
      coup_elems[1] = 16;
      coup_elems[2] = 19;
      coup_elems[3] = 14;
      coup_elems[4] = 16;
      coup_elems[5] = 14;
      coup_elems[6] = 19;
      coup_elems[7] = 21;
      HamiltHelper::prealloc_info(Delta,13,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(4,d_nnz[1]);
      EXPECT_EQ(4,o_nnz[1]);

      coup_elems[0] = 24;
      coup_elems[1] = 22;
      coup_elems[2] = 13;
      coup_elems[3] = 7;
      coup_elems[4] = 13;
      coup_elems[5] = 22;
      coup_elems[6] = 24;
      coup_elems[7] = 7;
      HamiltHelper::prealloc_info(Delta,14,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(6,o_nnz[2]);

      coup_elems[0] = 26;
      coup_elems[1] = 1;
      coup_elems[2] = 26;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,15,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[3]);
      EXPECT_EQ(4,o_nnz[3]);

      coup_elems[0] = 27;
      coup_elems[1] = 7;
      coup_elems[2] = 13;
      coup_elems[3] = 25;
      coup_elems[4] = 13;
      coup_elems[5] = 25;
      coup_elems[6] = 27;
      coup_elems[7] = 7;
      HamiltHelper::prealloc_info(Delta,16,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(6,o_nnz[4]);
      
      coup_elems[0] = 1;
      coup_elems[1] = 26;
      coup_elems[2] = 26;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,17,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[5]);
      EXPECT_EQ(4,o_nnz[5]);
    }

    if (mpi_rank == 3) {

      coup_elems[0] = 23;
      coup_elems[1] = 8;
      coup_elems[2] = 23;
      coup_elems[3] = 8;
      HamiltHelper::prealloc_info(Delta,18,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(4,o_nnz[0]);

      coup_elems[0] = 25;
      coup_elems[1] = 2;
      coup_elems[2] = 13;
      coup_elems[3] = 22;
      coup_elems[4] = 25;
      coup_elems[5] = 22;
      coup_elems[6] = 13;
      coup_elems[7] = 2;
      HamiltHelper::prealloc_info(Delta,19,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(6,o_nnz[1]);

      coup_elems[0] = 8;
      coup_elems[1] = 23;
      coup_elems[2] = 23;
      coup_elems[3] = 8;
      HamiltHelper::prealloc_info(Delta,20,gstart,gend,coup_elems,4,8,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[2]);
      EXPECT_EQ(4,o_nnz[2]);
      
      coup_elems[0] = 13;
      coup_elems[1] = 2;
      coup_elems[2] = 27;
      coup_elems[3] = 24;
      coup_elems[4] = 27;
      coup_elems[5] = 24;
      coup_elems[6] = 13;
      coup_elems[7] = 2;
      HamiltHelper::prealloc_info(Delta,21,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[3]);
      EXPECT_EQ(8,o_nnz[3]);

      coup_elems[0] = 3;
      coup_elems[1] = 14;
      coup_elems[2] = 19;
      coup_elems[3] = 10;
      coup_elems[4] = 19;
      coup_elems[5] = 14;
      coup_elems[6] = 10;
      coup_elems[7] = 3;
      HamiltHelper::prealloc_info(Delta,22,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[4]);
      EXPECT_EQ(6,o_nnz[4]);
    }

    if (mpi_rank == 4) {

      coup_elems[0] = 9;
      coup_elems[1] = 18;
      coup_elems[2] = 20;
      coup_elems[3] = 11;
      coup_elems[4] = 18;
      coup_elems[5] = 20;
      coup_elems[6] = 11;
      coup_elems[7] = 9;
      HamiltHelper::prealloc_info(Delta,23,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(8,o_nnz[0]);

      coup_elems[0] = 14;
      coup_elems[1] = 3;
      coup_elems[2] = 21;
      coup_elems[3] = 12;
      coup_elems[4] = 21;
      coup_elems[5] = 14;
      coup_elems[6] = 12;
      coup_elems[7] = 3;
      HamiltHelper::prealloc_info(Delta,24,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[1]);
      EXPECT_EQ(8,o_nnz[1]);

      coup_elems[0] = 10;
      coup_elems[1] = 19;
      coup_elems[2] = 5;
      coup_elems[3] = 16;
      coup_elems[4] = 19;
      coup_elems[5] = 16;
      coup_elems[6] = 10;
      coup_elems[7] = 5;
      HamiltHelper::prealloc_info(Delta,25,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[2]);
      EXPECT_EQ(8,o_nnz[2]);

      coup_elems[0] = 15;
      coup_elems[1] = 4;
      coup_elems[2] = 6;
      coup_elems[3] = 17;
      coup_elems[4] = 17;
      coup_elems[5] = 15;
      coup_elems[6] = 4;
      coup_elems[7] = 6;
      HamiltHelper::prealloc_info(Delta,26,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[3]);
      EXPECT_EQ(8,o_nnz[3]);

      coup_elems[0] = 16;
      coup_elems[1] = 5;
      coup_elems[2] = 12;
      coup_elems[3] = 21;
      coup_elems[4] = 21;
      coup_elems[5] = 16;
      coup_elems[6] = 12;
      coup_elems[7] = 5;
      HamiltHelper::prealloc_info(Delta,27,gstart,gend,coup_elems,8,0,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[4]);
      EXPECT_EQ(8,o_nnz[4]);
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

