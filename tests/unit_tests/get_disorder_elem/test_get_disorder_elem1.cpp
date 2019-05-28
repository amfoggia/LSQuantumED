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
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(-1.5, (double)HamiltHelper::get_disorder_elem(7,6,hi), 1e-13);
  EXPECT_NEAR(-1.166666666666666, (double)HamiltHelper::get_disorder_elem(11,6,hi), 1e-13);
  EXPECT_NEAR(-0.833333333333333, (double)HamiltHelper::get_disorder_elem(13,6,hi), 1e-13);
  EXPECT_NEAR(-0.5, (double)HamiltHelper::get_disorder_elem(14,6,hi), 1e-13);
  EXPECT_NEAR(-0.833333333333333, (double)HamiltHelper::get_disorder_elem(19,6,hi), 1e-13);
  EXPECT_NEAR(-0.5, (double)HamiltHelper::get_disorder_elem(21,6,hi), 1e-13);
  EXPECT_NEAR(-0.166666666666666, (double)HamiltHelper::get_disorder_elem(22,6,hi), 1e-13);
  EXPECT_NEAR(-0.166666666666666, (double)HamiltHelper::get_disorder_elem(25,6,hi), 1e-13);
  EXPECT_NEAR(0.166666666666666, (double)HamiltHelper::get_disorder_elem(26,6,hi), 1e-13);
  EXPECT_NEAR(0.5, (double)HamiltHelper::get_disorder_elem(28,6,hi), 1e-13);
  EXPECT_NEAR(-0.5, (double)HamiltHelper::get_disorder_elem(35,6,hi), 1e-13);
  EXPECT_NEAR(-0.166666666666666, (double)HamiltHelper::get_disorder_elem(37,6,hi), 1e-13);
  EXPECT_NEAR(0.166666666666666, (double)HamiltHelper::get_disorder_elem(38,6,hi), 1e-13);
  EXPECT_NEAR(0.166666666666666, (double)HamiltHelper::get_disorder_elem(41,6,hi), 1e-13);
  EXPECT_NEAR(0.5, (double)HamiltHelper::get_disorder_elem(42,6,hi), 1e-13);
  EXPECT_NEAR(0.833333333333333, (double)HamiltHelper::get_disorder_elem(44,6,hi), 1e-13);
  EXPECT_NEAR(0.5, (double)HamiltHelper::get_disorder_elem(49,6,hi), 1e-13);
  EXPECT_NEAR(0.833333333333333, (double)HamiltHelper::get_disorder_elem(50,6,hi), 1e-13);
  EXPECT_NEAR(1.166666666666666, (double)HamiltHelper::get_disorder_elem(52,6,hi), 1e-13);
  EXPECT_NEAR(1.5, (double)HamiltHelper::get_disorder_elem(56,6,hi), 1e-13);
}

TEST(GetElemDiag, chain1D61) {
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(-3.0/6.0, (double)HamiltHelper::get_disorder_elem(15,6,hi), 1e-13);
  EXPECT_NEAR(-1.0/6.0, (double)HamiltHelper::get_disorder_elem(23,6,hi), 1e-13);
  EXPECT_NEAR(1.0/6.0, (double)HamiltHelper::get_disorder_elem(27,6,hi), 1e-13);
  EXPECT_NEAR(3.0/6.0, (double)HamiltHelper::get_disorder_elem(29,6,hi), 1e-13);
  EXPECT_NEAR(5.0/6.0, (double)HamiltHelper::get_disorder_elem(30,6,hi), 1e-13);
  EXPECT_NEAR(1.0/6.0, (double)HamiltHelper::get_disorder_elem(39,6,hi), 1e-13);
  EXPECT_NEAR(3.0/6.0, (double)HamiltHelper::get_disorder_elem(43,6,hi), 1e-13);
  EXPECT_NEAR(5.0/6.0, (double)HamiltHelper::get_disorder_elem(45,6,hi), 1e-13);
  EXPECT_NEAR(7.0/6.0, (double)HamiltHelper::get_disorder_elem(46,6,hi), 1e-13);
  EXPECT_NEAR(5.0/6.0, (double)HamiltHelper::get_disorder_elem(51,6,hi), 1e-13);
  EXPECT_NEAR(7.0/6.0, (double)HamiltHelper::get_disorder_elem(53,6,hi), 1e-13);
  EXPECT_NEAR(9.0/6.0, (double)HamiltHelper::get_disorder_elem(54,6,hi), 1e-13);
  EXPECT_NEAR(9.0/6.0, (double)HamiltHelper::get_disorder_elem(57,6,hi), 1e-13);
  EXPECT_NEAR(11.0/6.0, (double)HamiltHelper::get_disorder_elem(58,6,hi), 1e-13);
  EXPECT_NEAR(13.0/6.0, (double)HamiltHelper::get_disorder_elem(60,6,hi), 1e-13);
}

TEST(GetElemDiag, chain1D6m1) {
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(-13.0/6.0, (double)HamiltHelper::get_disorder_elem(3,6,hi), 1e-13);
  EXPECT_NEAR(-11.0/6.0, (double)HamiltHelper::get_disorder_elem(5,6,hi), 1e-13);
  EXPECT_NEAR(-9.0/6.0, (double)HamiltHelper::get_disorder_elem(6,6,hi), 1e-13);
  EXPECT_NEAR(-9.0/6.0, (double)HamiltHelper::get_disorder_elem(9,6,hi), 1e-13);
  EXPECT_NEAR(-7.0/6.0, (double)HamiltHelper::get_disorder_elem(10,6,hi), 1e-13);
  EXPECT_NEAR(-5.0/6.0, (double)HamiltHelper::get_disorder_elem(12,6,hi), 1e-13);
  EXPECT_NEAR(-7.0/6.0, (double)HamiltHelper::get_disorder_elem(17,6,hi), 1e-13);
  EXPECT_NEAR(-5.0/6.0, (double)HamiltHelper::get_disorder_elem(18,6,hi), 1e-13);
  EXPECT_NEAR(-3.0/6.0, (double)HamiltHelper::get_disorder_elem(20,6,hi), 1e-13);
  EXPECT_NEAR(-1.0/6.0, (double)HamiltHelper::get_disorder_elem(24,6,hi), 1e-13);
  EXPECT_NEAR(-5.0/6.0, (double)HamiltHelper::get_disorder_elem(33,6,hi), 1e-13);
  EXPECT_NEAR(-3.0/6.0, (double)HamiltHelper::get_disorder_elem(34,6,hi), 1e-13);
  EXPECT_NEAR(-1.0/6.0, (double)HamiltHelper::get_disorder_elem(36,6,hi), 1e-13);
  EXPECT_NEAR(1.0/6.0, (double)HamiltHelper::get_disorder_elem(40,6,hi), 1e-13);
  EXPECT_NEAR(3.0/6.0, (double)HamiltHelper::get_disorder_elem(48,6,hi), 1e-13);
}

TEST(GetElemDiag, chain1D62) {
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(5.0/6.0, (double)HamiltHelper::get_disorder_elem(31,6,hi), 1e-13);
  EXPECT_NEAR(7.0/6.0, (double)HamiltHelper::get_disorder_elem(47,6,hi), 1e-13);
  EXPECT_NEAR(9.0/6.0, (double)HamiltHelper::get_disorder_elem(55,6,hi), 1e-13);
  EXPECT_NEAR(11.0/6.0, (double)HamiltHelper::get_disorder_elem(59,6,hi), 1e-13);
  EXPECT_NEAR(13.0/6.0, (double)HamiltHelper::get_disorder_elem(61,6,hi), 1e-13);
  EXPECT_NEAR(15.0/6.0, (double)HamiltHelper::get_disorder_elem(62,6,hi), 1e-13);
}

TEST(GetElemDiag, chain1D6m2) {
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(-15.0/6.0, (double)HamiltHelper::get_disorder_elem(1,6,hi), 1e-13);
  EXPECT_NEAR(-13.0/6.0, (double)HamiltHelper::get_disorder_elem(2,6,hi), 1e-13);
  EXPECT_NEAR(-11.0/6.0, (double)HamiltHelper::get_disorder_elem(4,6,hi), 1e-13);
  EXPECT_NEAR(-9.0/6.0, (double)HamiltHelper::get_disorder_elem(8,6,hi), 1e-13);
  EXPECT_NEAR(-7.0/6.0, (double)HamiltHelper::get_disorder_elem(16,6,hi), 1e-13);
  EXPECT_NEAR(-5.0/6.0, (double)HamiltHelper::get_disorder_elem(32,6,hi), 1e-13);
}

TEST(GetElemDiag, chain1D63) {
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(15.0/6.0, (double)HamiltHelper::get_disorder_elem(63,6,hi), 1e-13);
}

TEST(GetElemDiag, chain1D6m3) {
  PetscReal * hi;
  PetscCalloc1(env.nspins, &hi);
  for (PetscInt i = 0; i < env.nspins; ++i)
    hi[i] = PetscReal(i)/env.nspins; 

  EXPECT_NEAR(-15.0/6.0, (double)HamiltHelper::get_disorder_elem(0,6,hi), 1e-13);
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

