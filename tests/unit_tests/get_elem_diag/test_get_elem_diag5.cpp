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

TEST(GetElemDiag, square2D125) {
  honeycomb2D l1{env,2,4};
  
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,63,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,95,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,111,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,119,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,123,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,125,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,126,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,159,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,175,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,183,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,187,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,189,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,190,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,207,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,215,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,219,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,221,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,222,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,231,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,235,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,237,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,238,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,243,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,245,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,246,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,249,8,l1));
  EXPECT_EQ(0.0,HamiltHelper::get_elem_diag(neigh_type::nn,250,8,l1));
  EXPECT_EQ(4.0,HamiltHelper::get_elem_diag(neigh_type::nn,252,8,l1));
  
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

