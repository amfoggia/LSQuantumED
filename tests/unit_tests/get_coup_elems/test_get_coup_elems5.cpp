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

TEST(GetCoupElems, honeycomb2D125NN) {
  Basis b{env,2};
  honeycomb2D lat{env,2,4};
  PetscInt ncol;
  PetscErrorCode ierr = 0;
  PetscInt * coup_elems;
  PetscCalloc1(b.nspins-1, &coup_elems);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,0,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(4,coup_elems[0]);
  EXPECT_EQ(1,coup_elems[1]);
  EXPECT_EQ(8,coup_elems[2]);
  EXPECT_EQ(9,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,1,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(2,coup_elems[0]);
  EXPECT_EQ(0,coup_elems[1]);
  EXPECT_EQ(5,coup_elems[2]);
  EXPECT_EQ(7,coup_elems[3]);
  EXPECT_EQ(13,coup_elems[4]);
  EXPECT_EQ(14,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,2,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(6,coup_elems[0]);
  EXPECT_EQ(1,coup_elems[1]);
  EXPECT_EQ(8,coup_elems[2]);
  EXPECT_EQ(18,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,3,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(4,coup_elems[0]);
  EXPECT_EQ(6,coup_elems[1]);
  EXPECT_EQ(9,coup_elems[2]);
  EXPECT_EQ(18,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,4,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(5,coup_elems[0]);
  EXPECT_EQ(3,coup_elems[1]);
  EXPECT_EQ(0,coup_elems[2]);
  EXPECT_EQ(10,coup_elems[3]);
  EXPECT_EQ(19,coup_elems[4]);
  EXPECT_EQ(22,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,5,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(6,coup_elems[0]);
  EXPECT_EQ(4,coup_elems[1]);
  EXPECT_EQ(1,coup_elems[2]);
  EXPECT_EQ(11,coup_elems[3]);
  EXPECT_EQ(20,coup_elems[4]);
  EXPECT_EQ(23,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,6,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(5,coup_elems[0]);
  EXPECT_EQ(2,coup_elems[1]);
  EXPECT_EQ(3,coup_elems[2]);
  EXPECT_EQ(12,coup_elems[3]);
  EXPECT_EQ(21,coup_elems[4]);
  EXPECT_EQ(24,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,7,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(15,coup_elems[0]);
  EXPECT_EQ(8,coup_elems[1]);
  EXPECT_EQ(11,coup_elems[2]);
  EXPECT_EQ(1,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,8,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(12,coup_elems[0]);
  EXPECT_EQ(19,coup_elems[1]);
  EXPECT_EQ(7,coup_elems[2]);
  EXPECT_EQ(13,coup_elems[3]);
  EXPECT_EQ(2,coup_elems[4]);
  EXPECT_EQ(0,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,9,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(10,coup_elems[0]);
  EXPECT_EQ(22,coup_elems[1]);
  EXPECT_EQ(12,coup_elems[2]);
  EXPECT_EQ(14,coup_elems[3]);
  EXPECT_EQ(3,coup_elems[4]);
  EXPECT_EQ(0,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,10,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(11,coup_elems[0]);
  EXPECT_EQ(9,coup_elems[1]);
  EXPECT_EQ(15,coup_elems[2]);
  EXPECT_EQ(4,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,11,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(12,coup_elems[0]);
  EXPECT_EQ(10,coup_elems[1]);
  EXPECT_EQ(25,coup_elems[2]);
  EXPECT_EQ(16,coup_elems[3]);
  EXPECT_EQ(7,coup_elems[4]);
  EXPECT_EQ(5,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,12,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(11,coup_elems[0]);
  EXPECT_EQ(8,coup_elems[1]);
  EXPECT_EQ(26,coup_elems[2]);
  EXPECT_EQ(9,coup_elems[3]);
  EXPECT_EQ(17,coup_elems[4]);
  EXPECT_EQ(6,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,13,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(17,coup_elems[0]);
  EXPECT_EQ(8,coup_elems[1]);
  EXPECT_EQ(20,coup_elems[2]);
  EXPECT_EQ(1,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,14,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(15,coup_elems[0]);
  EXPECT_EQ(17,coup_elems[1]);
  EXPECT_EQ(18,coup_elems[2]);
  EXPECT_EQ(9,coup_elems[3]);
  EXPECT_EQ(23,coup_elems[4]);
  EXPECT_EQ(1,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,15,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(16,coup_elems[0]);
  EXPECT_EQ(14,coup_elems[1]);
  EXPECT_EQ(7,coup_elems[2]);
  EXPECT_EQ(19,coup_elems[3]);
  EXPECT_EQ(10,coup_elems[4]);
  EXPECT_EQ(25,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,16,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(17,coup_elems[0]);
  EXPECT_EQ(15,coup_elems[1]);
  EXPECT_EQ(20,coup_elems[2]);
  EXPECT_EQ(11,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,17,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(16,coup_elems[0]);
  EXPECT_EQ(13,coup_elems[1]);
  EXPECT_EQ(14,coup_elems[2]);
  EXPECT_EQ(21,coup_elems[3]);
  EXPECT_EQ(12,coup_elems[4]);
  EXPECT_EQ(27,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,18,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(24,coup_elems[0]);
  EXPECT_EQ(19,coup_elems[1]);
  EXPECT_EQ(21,coup_elems[2]);
  EXPECT_EQ(14,coup_elems[3]);
  EXPECT_EQ(3,coup_elems[4]);
  EXPECT_EQ(2,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,19,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(26,coup_elems[0]);
  EXPECT_EQ(20,coup_elems[1]);
  EXPECT_EQ(18,coup_elems[2]);
  EXPECT_EQ(8,coup_elems[3]);
  EXPECT_EQ(15,coup_elems[4]);
  EXPECT_EQ(4,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,20,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(21,coup_elems[0]);
  EXPECT_EQ(27,coup_elems[1]);
  EXPECT_EQ(19,coup_elems[2]);
  EXPECT_EQ(16,coup_elems[3]);
  EXPECT_EQ(13,coup_elems[4]);
  EXPECT_EQ(5,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,21,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(20,coup_elems[0]);
  EXPECT_EQ(18,coup_elems[1]);
  EXPECT_EQ(17,coup_elems[2]);
  EXPECT_EQ(6,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,22,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(23,coup_elems[0]);
  EXPECT_EQ(9,coup_elems[1]);
  EXPECT_EQ(26,coup_elems[2]);
  EXPECT_EQ(4,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,23,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(24,coup_elems[0]);
  EXPECT_EQ(22,coup_elems[1]);
  EXPECT_EQ(25,coup_elems[2]);
  EXPECT_EQ(27,coup_elems[3]);
  EXPECT_EQ(14,coup_elems[4]);
  EXPECT_EQ(5,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,24,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(23,coup_elems[0]);
  EXPECT_EQ(18,coup_elems[1]);
  EXPECT_EQ(26,coup_elems[2]);
  EXPECT_EQ(6,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,25,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(26,coup_elems[0]);
  EXPECT_EQ(23,coup_elems[1]);
  EXPECT_EQ(11,coup_elems[2]);
  EXPECT_EQ(15,coup_elems[3]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,26,lat,ncol,coup_elems);
  ASSERT_EQ(6,ncol);
  EXPECT_EQ(25,coup_elems[0]);
  EXPECT_EQ(19,coup_elems[1]);
  EXPECT_EQ(27,coup_elems[2]);
  EXPECT_EQ(24,coup_elems[3]);
  EXPECT_EQ(12,coup_elems[4]);
  EXPECT_EQ(22,coup_elems[5]);

  HamiltHelper::get_coup_elems(neigh_type::nn,b,27,lat,ncol,coup_elems);
  ASSERT_EQ(4,ncol);
  EXPECT_EQ(20,coup_elems[0]);
  EXPECT_EQ(26,coup_elems[1]);
  EXPECT_EQ(23,coup_elems[2]);
  EXPECT_EQ(17,coup_elems[3]);

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
