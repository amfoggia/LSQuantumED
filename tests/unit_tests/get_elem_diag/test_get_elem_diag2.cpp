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

TEST(GetElemDiag, square2D82NN) {
  square2D l1{env,2,4};
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,63,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,95,8,l1));
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,111,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,119,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,123,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,125,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,126,8,l1));
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,159,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,175,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,183,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,187,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,189,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,190,8,l1));
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,207,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,215,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,219,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,221,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,222,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,231,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,235,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,237,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,238,8,l1));
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,243,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,245,8,l1));
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,246,8,l1));
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,249,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nn,250,8,l1));
  EXPECT_EQ(4,HamiltHelper::get_elem_diag(neigh_type::nn,252,8,l1));
}

TEST(GetElemDiag, square2D82NNN) {
  square2D l1{env,2,4};
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,63,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,95,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,111,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,119,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,123,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,125,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,126,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,159,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,175,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,183,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,187,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,189,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,190,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,207,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,215,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,219,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,221,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,222,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,231,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,235,8,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,237,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,238,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,243,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,245,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,246,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,249,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,250,8,l1));
  EXPECT_EQ(0,HamiltHelper::get_elem_diag(neigh_type::nnn,252,8,l1));
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

