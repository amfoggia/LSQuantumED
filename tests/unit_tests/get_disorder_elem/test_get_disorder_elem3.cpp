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
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(44.0/12.0,(double)HamiltHelper::get_disorder_elem(2047,12,hi), 1e-13);
  EXPECT_NEAR(46.0/12.0,(double)HamiltHelper::get_disorder_elem(3071,12,hi), 1e-13);
  EXPECT_NEAR(48.0/12.0,(double)HamiltHelper::get_disorder_elem(3583,12,hi), 1e-13);
  EXPECT_NEAR(50.0/12.0,(double)HamiltHelper::get_disorder_elem(3839,12,hi), 1e-13);
  EXPECT_NEAR(52.0/12.0,(double)HamiltHelper::get_disorder_elem(3967,12,hi), 1e-13);
  EXPECT_NEAR(54.0/12.0,(double)HamiltHelper::get_disorder_elem(4031,12,hi), 1e-13);
  EXPECT_NEAR(56.0/12.0,(double)HamiltHelper::get_disorder_elem(4063,12,hi), 1e-13);
  EXPECT_NEAR(58.0/12.0,(double)HamiltHelper::get_disorder_elem(4079,12,hi), 1e-13);
  EXPECT_NEAR(60.0/12.0,(double)HamiltHelper::get_disorder_elem(4087,12,hi), 1e-13);
  EXPECT_NEAR(62.0/12.0,(double)HamiltHelper::get_disorder_elem(4091,12,hi), 1e-13);
  EXPECT_NEAR(64.0/12.0,(double)HamiltHelper::get_disorder_elem(4093,12,hi), 1e-13);
  EXPECT_NEAR(66.0/12.0,(double)HamiltHelper::get_disorder_elem(4094,12,hi), 1e-13);
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

