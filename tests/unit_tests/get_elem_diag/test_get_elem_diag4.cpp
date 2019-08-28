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

TEST(GetElemDiag, square2D125) {
  honeycomb2D l1{env,6,2};
  
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,2047,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,3071,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,3583,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,3839,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,3967,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,4031,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,4063,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,4079,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,4087,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,4091,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,4093,12,l1));
  EXPECT_EQ(12.0,HamiltHelper::get_elem_diag(neigh_type::nn,4094,12,l1));

  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,2047,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3071,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3583,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3839,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3967,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4031,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4063,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4079,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4087,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4091,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4093,12,l1));
  // EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4094,12,l1));
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

