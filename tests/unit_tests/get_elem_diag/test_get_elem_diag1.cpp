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


TEST(GetElemDiag, chain1D) {
  chain1D l1{env};
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,7,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,11,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,13,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,14,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,19,6,l1));
  EXPECT_EQ(-6, HamiltHelper::get_elem_diag(neigh_type::nn,21,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,22,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,25,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,26,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,28,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,35,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,37,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,38,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,41,6,l1));
  EXPECT_EQ(-6, HamiltHelper::get_elem_diag(neigh_type::nn,42,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,44,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,49,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,50,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,52,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,56,6,l1));

  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,7,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,11,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,13,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,14,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,19,6,l1));
  EXPECT_EQ(6, HamiltHelper::get_elem_diag(neigh_type::nnn,21,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,22,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,25,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,26,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,28,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,35,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,37,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,38,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,41,6,l1));
  EXPECT_EQ(6, HamiltHelper::get_elem_diag(neigh_type::nnn,42,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,44,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,49,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,50,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,52,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,56,6,l1));
}

TEST(GetElemDiag, chain1D61) {
  chain1D l1{env};
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,15,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,23,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,27,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,29,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,30,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,39,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,43,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,45,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,46,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,51,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,53,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,54,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,57,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,58,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,60,6,l1));

  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,15,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,23,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,27,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,29,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,30,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,39,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,43,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,45,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,46,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,51,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,53,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,54,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,57,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,58,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,60,6,l1));
}

TEST(GetElemDiag, chain1D6m1) {
  chain1D l1{env};
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,3,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,5,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,6,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,9,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,10,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,12,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,17,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,18,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,20,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,24,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,33,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,34,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,36,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nn,40,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,48,6,l1));

  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,3,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,5,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,6,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,9,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,10,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,12,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,17,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,18,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,20,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,24,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,33,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,34,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,36,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,40,6,l1));
  EXPECT_EQ(-2, HamiltHelper::get_elem_diag(neigh_type::nnn,48,6,l1));
}

TEST(GetElemDiag, chain1D62) {
  chain1D l1{env};
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,31,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,47,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,55,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,59,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,61,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,62,6,l1));

  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,31,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,47,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,55,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,59,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,61,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,62,6,l1));
}

TEST(GetElemDiag, chain1D6m2) {
  chain1D l1{env};
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,1,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,2,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,4,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,8,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,16,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nn,32,6,l1));

  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,1,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,2,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,4,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,8,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,16,6,l1));
  EXPECT_EQ(2, HamiltHelper::get_elem_diag(neigh_type::nnn,32,6,l1));
}

TEST(GetElemDiag, chain1D63) {
  chain1D l1{env};
  EXPECT_EQ(6, HamiltHelper::get_elem_diag(neigh_type::nn,63,6,l1));
}

TEST(GetElemDiag, chain1D6m3) {
  chain1D l1{env};
  EXPECT_EQ(6, HamiltHelper::get_elem_diag(neigh_type::nn,0,6,l1));
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

