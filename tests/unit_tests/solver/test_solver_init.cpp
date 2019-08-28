#include <gtest/gtest.h>
#include "hamiltonian.hpp"
#include "solver.hpp"

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Solver namespace\n\n";
int _argc;
char ** _argv;

Environment env{_argc,_argv,help};

class HamiltonianTestEnv : public ::testing::Environment {
protected:
  virtual void TearDown() {
    ::testing::internal::CaptureStdout();
    EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};


TEST(HamiltonianTest, SolverInit) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  env.nspins = 6;
  Basis b{env,0};
  chain1D lat{env};
  Hamiltonian<chain1D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;

  EPS solver;
  EPSType type;
  EPSProblemType ptype;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);
  
  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,1));

  ierr = EPSGetType(solver,&type); ASSERT_EQ(0,ierr);
  EXPECT_STREQ("krylovschur", type);
  
  ierr = EPSGetProblemType(solver,&ptype); ASSERT_EQ(0,ierr);
  EXPECT_EQ(1, ptype);

  EXPECT_EQ(0,Solver::SolverClean(solver));
}


int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&_argc, _argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}
