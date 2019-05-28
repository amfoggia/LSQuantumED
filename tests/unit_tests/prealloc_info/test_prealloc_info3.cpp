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

Environment env{_argc,_argv,12,help};

class HamiltonianTestEnv : public ::testing::Environment {
protected:
  virtual void TearDown() {
    ::testing::internal::CaptureStdout();
    EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};

TEST(PreallocInfo, square2D125) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 1.0;

    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 12;
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

      coup_elems[0] = 4;
      coup_elems[1] = 1;
      coup_elems[2] = 3;
      coup_elems[3] = 8;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[0]);
      EXPECT_EQ(2,o_nnz[0]);

      coup_elems[0] = 5;
      coup_elems[1] = 2;
      coup_elems[2] = 0;
      coup_elems[3] = 9;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 6;
      coup_elems[1] = 3;
      coup_elems[2] = 1;
      coup_elems[3] = 10;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(2,o_nnz[2]);

      coup_elems[0] = 7;
      coup_elems[1] = 2;
      coup_elems[2] = 11;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);    
    }
    
    if (mpi_rank == 1) {

      coup_elems[0] = 8;
      coup_elems[1] = 5;
      coup_elems[2] = 7;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[0]);
      EXPECT_EQ(2,o_nnz[0]);

      coup_elems[0] = 9;
      coup_elems[1] = 6;
      coup_elems[2] = 4;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 10;
      coup_elems[1] = 7;
      coup_elems[2] = 5;
      coup_elems[3] = 2;
      HamiltHelper::prealloc_info(Delta,6,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(2,o_nnz[2]);

      coup_elems[0] = 11;
      coup_elems[1] = 6;
      coup_elems[2] = 3;
      coup_elems[3] = 4;
      HamiltHelper::prealloc_info(Delta,7,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 9;
      coup_elems[1] = 11;
      coup_elems[2] = 4;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,8,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[0]);
      EXPECT_EQ(2,o_nnz[0]);

      coup_elems[0] = 10;
      coup_elems[1] = 8;
      coup_elems[2] = 5;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,9,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[1]);
      EXPECT_EQ(2,o_nnz[1]);

      coup_elems[0] = 11;
      coup_elems[1] = 9;
      coup_elems[2] = 6;
      coup_elems[3] = 2;
      HamiltHelper::prealloc_info(Delta,10,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[2]);
      EXPECT_EQ(2,o_nnz[2]);

      coup_elems[0] = 10;
      coup_elems[1] = 7;
      coup_elems[2] = 8;
      coup_elems[3] = 3;
      HamiltHelper::prealloc_info(Delta,11,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(2,d_nnz[3]);
      EXPECT_EQ(2,o_nnz[3]);
    }

    ierr = PetscFree(coup_elems); EXPECT_EQ(0,ierr);
    ierr = PetscFree(d_nnz); EXPECT_EQ(0,ierr);
    ierr = PetscFree(o_nnz); EXPECT_EQ(0,ierr);
  }
}

TEST(PreallocInfo, square2D125_J2) {
  {
    PetscErrorCode ierr = 0;

    PetscReal Delta = 0.0;

    PetscMPIInt mpi_rank, mpi_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
  
    PetscInt size = 12;
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

      coup_elems[0] = 7;
      coup_elems[1] = 5;
      coup_elems[2] = 11;
      coup_elems[3] = 9;
      HamiltHelper::prealloc_info(Delta,0,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(4,o_nnz[0]);

      coup_elems[0] = 6;
      coup_elems[1] = 4;
      coup_elems[2] = 8;
      coup_elems[3] = 10;
      HamiltHelper::prealloc_info(Delta,1,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[1]);
      EXPECT_EQ(4,o_nnz[1]);

      coup_elems[0] = 7;
      coup_elems[1] = 5;
      coup_elems[2] = 9;
      coup_elems[3] = 11;
      HamiltHelper::prealloc_info(Delta,2,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[2]);
      EXPECT_EQ(4,o_nnz[2]);

      coup_elems[0] = 6;
      coup_elems[1] = 4;
      coup_elems[2] = 10;
      coup_elems[3] = 8;
      HamiltHelper::prealloc_info(Delta,3,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[3]);
      EXPECT_EQ(4,o_nnz[3]);    
    }
    
    if (mpi_rank == 1) {

      coup_elems[0] = 11;
      coup_elems[1] = 9;
      coup_elems[2] = 3;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,4,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(4,o_nnz[0]);

      coup_elems[0] = 10;
      coup_elems[1] = 8;
      coup_elems[2] = 0;
      coup_elems[3] = 2;
      HamiltHelper::prealloc_info(Delta,5,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[1]);
      EXPECT_EQ(4,o_nnz[1]);

      coup_elems[0] = 11;
      coup_elems[1] = 9;
      coup_elems[2] = 1;
      coup_elems[3] = 3;
      HamiltHelper::prealloc_info(Delta,6,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[2]);
      EXPECT_EQ(4,o_nnz[2]);

      coup_elems[0] = 10;
      coup_elems[1] = 8;
      coup_elems[2] = 2;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,7,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[3]);
      EXPECT_EQ(4,o_nnz[3]);
    }

    if (mpi_rank == 2) {

      coup_elems[0] = 7;
      coup_elems[1] = 5;
      coup_elems[2] = 3;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,8,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[0]);
      EXPECT_EQ(4,o_nnz[0]);

      coup_elems[0] = 4;
      coup_elems[1] = 6;
      coup_elems[2] = 2;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,9,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[1]);
      EXPECT_EQ(4,o_nnz[1]);

      coup_elems[0] = 5;
      coup_elems[1] = 7;
      coup_elems[2] = 3;
      coup_elems[3] = 1;
      HamiltHelper::prealloc_info(Delta,10,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[2]);
      EXPECT_EQ(4,o_nnz[2]);

      coup_elems[0] = 6;
      coup_elems[1] = 4;
      coup_elems[2] = 2;
      coup_elems[3] = 0;
      HamiltHelper::prealloc_info(Delta,11,gstart,gend,coup_elems,4,16,d_nnz,o_nnz);
      EXPECT_EQ(0,d_nnz[3]);
      EXPECT_EQ(4,o_nnz[3]);
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

