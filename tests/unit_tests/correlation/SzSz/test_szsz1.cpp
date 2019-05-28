#include "hamiltonian.hpp"
#include "solver.hpp"
#include "correlation.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Correlation-SzSz class\n\n";
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


TEST(SzSzCorrelation, chain1D) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  Basis b{env,0};
  chain1D lat{env};
  Hamiltonian<chain1D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;
  SzSz corr{env,b};
  PetscScalar * exp_val;
  PetscCalloc1(b.nspins-1, &exp_val);

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

#ifdef DISORDER
  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
#endif
  
  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,1));
 
  ierr = Solver::solve_lanczos(env,solver,nconv); EXPECT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); EXPECT_EQ(0,ierr);

  ASSERT_GE(nconv,1);
  if (nconv >= 1) {
    ierr = Solver::solution(solver, 0, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-3.6510934089371827, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  //  ierr = VecView(xr, PETSC_VIEWER_STDOUT_WORLD);

  for (PetscInt i = 0; i < b.nspins-1; ++i)
    exp_val[i] = corr.ExpectVal(env,0,i+1,xr);
  
  EXPECT_NEAR(-0.15212889203236324, PetscRealPart(exp_val[0]), 5e-10);
  EXPECT_NEAR(0.0652593013325114, PetscRealPart(exp_val[1]), 5e-10);
  EXPECT_NEAR(-0.0629842737332784, PetscRealPart(exp_val[2]), 5e-10);
  EXPECT_NEAR(0.0497077289039414, PetscRealPart(exp_val[3]), 5e-10);
  EXPECT_NEAR(-0.0629842737728710, PetscRealPart(exp_val[4]), 5e-10);
  EXPECT_NEAR(0.0652593013587580, PetscRealPart(exp_val[5]), 5e-10);
  EXPECT_NEAR(-0.1521288920566983, PetscRealPart(exp_val[6]), 5e-10);

  // for (PetscInt i = 0; i < b.nspins-1; ++i)
  //   PetscPrintf(PETSC_COMM_WORLD, "exp_val[%d]: %.16f\n", i, exp_val[i]);
  
  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0, Solver::SolverClean(solver));
  EXPECT_EQ(0, VecDestroy(&xr));
  EXPECT_EQ(0, PetscFree(exp_val));
  
}


int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

