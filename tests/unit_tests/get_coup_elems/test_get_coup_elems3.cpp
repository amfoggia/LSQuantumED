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

TEST(GetCoupElems, square2D125NN) {
  Basis b{env,5};
  square2D lat{env,3,4};
  PetscInt ncol;
  PetscErrorCode ierr = 0;
  PetscInt * coup_elems;
  PetscCalloc1(b.nspins-1, &coup_elems);

  if (b.mpi_rank == 0) {
    HamiltHelper::get_coup_elems(neigh_type::nn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(4,coup_elems[0]);
    EXPECT_EQ(1,coup_elems[1]);
    EXPECT_EQ(3,coup_elems[2]);
    EXPECT_EQ(8,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(5,coup_elems[0]);
    EXPECT_EQ(2,coup_elems[1]);
    EXPECT_EQ(0,coup_elems[2]);
    EXPECT_EQ(9,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nn,b,2,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(6,coup_elems[0]);
    EXPECT_EQ(3,coup_elems[1]);
    EXPECT_EQ(1,coup_elems[2]);
    EXPECT_EQ(10,coup_elems[3]);
  }

  else if (b.mpi_rank == 1) {
    HamiltHelper::get_coup_elems(neigh_type::nn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(7,coup_elems[0]);
    EXPECT_EQ(2,coup_elems[1]);
    EXPECT_EQ(11,coup_elems[2]);
    EXPECT_EQ(0,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(8,coup_elems[0]);
    EXPECT_EQ(5,coup_elems[1]);
    EXPECT_EQ(7,coup_elems[2]);
    EXPECT_EQ(0,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nn,b,2,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(9,coup_elems[0]);
    EXPECT_EQ(6,coup_elems[1]);
    EXPECT_EQ(4,coup_elems[2]);
    EXPECT_EQ(1,coup_elems[3]);
  }

  else if (b.mpi_rank == 2) {
    HamiltHelper::get_coup_elems(neigh_type::nn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(10,coup_elems[0]);
    EXPECT_EQ(7,coup_elems[1]);
    EXPECT_EQ(5,coup_elems[2]);
    EXPECT_EQ(2,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(11,coup_elems[0]);
    EXPECT_EQ(6,coup_elems[1]);
    EXPECT_EQ(3,coup_elems[2]);
    EXPECT_EQ(4,coup_elems[3]);
  }

  else if (b.mpi_rank == 3) {
    HamiltHelper::get_coup_elems(neigh_type::nn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(9,coup_elems[0]);
    EXPECT_EQ(11,coup_elems[1]);
    EXPECT_EQ(4,coup_elems[2]);
    EXPECT_EQ(0,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(10,coup_elems[0]);
    EXPECT_EQ(8,coup_elems[1]);
    EXPECT_EQ(5,coup_elems[2]);
    EXPECT_EQ(1,coup_elems[3]);
  }

  else {
    HamiltHelper::get_coup_elems(neigh_type::nn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(11,coup_elems[0]);
    EXPECT_EQ(9,coup_elems[1]);
    EXPECT_EQ(6,coup_elems[2]);
    EXPECT_EQ(2,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(10,coup_elems[0]);
    EXPECT_EQ(7,coup_elems[1]);
    EXPECT_EQ(8,coup_elems[2]);
    EXPECT_EQ(3,coup_elems[3]);
  }
  
  ierr = PetscFree(coup_elems);
  EXPECT_EQ(0,ierr);
}

TEST(GetCoupElems, square2D125NNN) {
  Basis b{env,5};
  square2D lat{env,3,4};
  PetscInt ncol;
  PetscErrorCode ierr = 0;
  PetscInt * coup_elems;
  PetscCalloc1(b.nspins-1, &coup_elems);

  if (b.mpi_rank == 0) {
    HamiltHelper::get_coup_elems(neigh_type::nnn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(7,coup_elems[0]);
    EXPECT_EQ(5,coup_elems[1]);
    EXPECT_EQ(11,coup_elems[2]);
    EXPECT_EQ(9,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nnn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(6,coup_elems[0]);
    EXPECT_EQ(4,coup_elems[1]);
    EXPECT_EQ(8,coup_elems[2]);
    EXPECT_EQ(10,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nnn,b,2,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(7,coup_elems[0]);
    EXPECT_EQ(5,coup_elems[1]);
    EXPECT_EQ(9,coup_elems[2]);
    EXPECT_EQ(11,coup_elems[3]);
  }

  else if (b.mpi_rank == 1) {
    HamiltHelper::get_coup_elems(neigh_type::nnn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(6,coup_elems[0]);
    EXPECT_EQ(4,coup_elems[1]);
    EXPECT_EQ(10,coup_elems[2]);
    EXPECT_EQ(8,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nnn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(11,coup_elems[0]);
    EXPECT_EQ(9,coup_elems[1]);
    EXPECT_EQ(3,coup_elems[2]);
    EXPECT_EQ(1,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nnn,b,2,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(10,coup_elems[0]);
    EXPECT_EQ(8,coup_elems[1]);
    EXPECT_EQ(0,coup_elems[2]);
    EXPECT_EQ(2,coup_elems[3]);
  }

  else if (b.mpi_rank == 2) {
    HamiltHelper::get_coup_elems(neigh_type::nnn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(11,coup_elems[0]);
    EXPECT_EQ(9,coup_elems[1]);
    EXPECT_EQ(1,coup_elems[2]);
    EXPECT_EQ(3,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nnn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(10,coup_elems[0]);
    EXPECT_EQ(8,coup_elems[1]);
    EXPECT_EQ(2,coup_elems[2]);
    EXPECT_EQ(0,coup_elems[3]);
  }

  else if (b.mpi_rank == 3) {
    HamiltHelper::get_coup_elems(neigh_type::nnn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(7,coup_elems[0]);
    EXPECT_EQ(5,coup_elems[1]);
    EXPECT_EQ(3,coup_elems[2]);
    EXPECT_EQ(1,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nnn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(4,coup_elems[0]);
    EXPECT_EQ(6,coup_elems[1]);
    EXPECT_EQ(2,coup_elems[2]);
    EXPECT_EQ(0,coup_elems[3]);
  }

  else {
    HamiltHelper::get_coup_elems(neigh_type::nnn,b,0,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(5,coup_elems[0]);
    EXPECT_EQ(7,coup_elems[1]);
    EXPECT_EQ(3,coup_elems[2]);
    EXPECT_EQ(1,coup_elems[3]);

    HamiltHelper::get_coup_elems(neigh_type::nnn,b,1,lat,ncol,coup_elems);
    ASSERT_EQ(4,ncol);
    EXPECT_EQ(6,coup_elems[0]);
    EXPECT_EQ(4,coup_elems[1]);
    EXPECT_EQ(2,coup_elems[2]);
    EXPECT_EQ(0,coup_elems[3]);
  }
  
  ierr = PetscFree(coup_elems);
  EXPECT_EQ(0,ierr);
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}
