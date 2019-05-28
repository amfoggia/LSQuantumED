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
  square2D l1{env,3,4};
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,2047,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,3071,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,3583,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,3839,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,3967,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,4031,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,4063,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,4079,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,4087,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,4091,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,4093,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nn,4094,12,l1));

  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,2047,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3071,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3583,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3839,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,3967,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4031,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4063,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4079,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4087,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4091,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4093,12,l1));
  EXPECT_EQ(16,HamiltHelper::get_elem_diag(neigh_type::nnn,4094,12,l1));
}

TEST(GetElemDiag, square2D124NN) {
  square2D l1{env,3,4};
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,1023,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,1535,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,1791,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,1919,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,1983,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,2015,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,2031,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,2039,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,2043,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,2045,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,2046,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,2559,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,2815,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,2943,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3007,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3039,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3055,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3063,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3067,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3069,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3070,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3327,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3455,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3519,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3551,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3567,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3575,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3579,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3581,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3582,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3711,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3775,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3807,12,l1));
  
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3823,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3831,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3835,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3837,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3838,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3903,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3935,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3951,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3959,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3963,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3965,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,3966,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,3999,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4015,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4023,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4027,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4029,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4030,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4047,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4055,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4059,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4061,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4062,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4071,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4075,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4077,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4078,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4083,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4085,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4086,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4089,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nn,4090,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nn,4092,12,l1));
}

TEST(GetElemDiag, square2D124NNN) {
  square2D l1{env,3,4};
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,1023,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,1535,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,1791,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,1919,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,1983,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,2015,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,2031,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,2039,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,2043,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,2045,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,2046,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,2559,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,2815,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,2943,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3007,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3039,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3055,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3063,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3067,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3069,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3070,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3327,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3455,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3519,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3551,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3567,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3575,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3579,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3581,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3582,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3711,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3775,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3807,12,l1));
  
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3823,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3831,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3835,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3837,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3838,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3903,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3935,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3951,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3959,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3963,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3965,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,3966,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,3999,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4015,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,4023,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4027,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,4029,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4030,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4047,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4055,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,4059,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4061,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,4062,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,4071,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4075,12,l1));
  EXPECT_EQ(12,HamiltHelper::get_elem_diag(neigh_type::nnn,4077,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4078,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4083,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4085,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4086,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4089,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4090,12,l1));
  EXPECT_EQ(8,HamiltHelper::get_elem_diag(neigh_type::nnn,4092,12,l1));
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

