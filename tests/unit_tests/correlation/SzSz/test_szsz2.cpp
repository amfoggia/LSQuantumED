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

Environment env{_argc,_argv,16,help};

class HamiltonianTestEnv : public ::testing::Environment {
protected:

  virtual void TearDown() {
    ::testing::internal::CaptureStdout();
    EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};


TEST(SzSzCorrelation, square2D_16spins) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  Basis b{env,0};
  square2D lat{env,4,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;
  SzSz corr{env,b};
  PetscScalar * exp_val;
  PetscCalloc1(b.nspins, &exp_val);

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);
  
  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,1));
 
  ierr = Solver::solve(env,solver,nconv); EXPECT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); EXPECT_EQ(0,ierr);

  ASSERT_GE(nconv,1);
  if (nconv >= 1) {
    ierr = Solver::solution(solver, 0, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-11.228483208428859, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  //  ierr = VecView(xr, PETSC_VIEWER_STDOUT_WORLD);

  for (PetscInt i = 0; i < b.nspins; ++i)
    exp_val[i] = corr.ExpectVal(env,0,i,xr);
  
  EXPECT_NEAR(0.2500000000000000, PetscRealPart(exp_val[0]), 5e-10);
  EXPECT_NEAR(-0.1169633667407941, PetscRealPart(exp_val[1]), 5e-10);
  EXPECT_NEAR(0.0712550951102676, PetscRealPart(exp_val[2]), 5e-10);
  EXPECT_NEAR(-0.1169633667170021, PetscRealPart(exp_val[3]), 5e-10);
  EXPECT_NEAR(-0.1169633667322925, PetscRealPart(exp_val[4]), 5e-10);
  EXPECT_NEAR(0.0712550951091117, PetscRealPart(exp_val[5]), 5e-10);
  EXPECT_NEAR(-0.0673880572985575, PetscRealPart(exp_val[6]), 5e-10);
  EXPECT_NEAR(0.0712550951125746, PetscRealPart(exp_val[7]), 5e-10);
  EXPECT_NEAR(0.0712550951162089, PetscRealPart(exp_val[8]), 5e-10);
  EXPECT_NEAR(-0.0673880573063124, PetscRealPart(exp_val[9]), 5e-10);
  EXPECT_NEAR(0.0598751254796498, PetscRealPart(exp_val[10]), 5e-10);
  EXPECT_NEAR(-0.0673880572996231, PetscRealPart(exp_val[11]), 5e-10);
  EXPECT_NEAR(-0.1169633667543916, PetscRealPart(exp_val[12]), 5e-10);
  EXPECT_NEAR(0.0712550951182559, PetscRealPart(exp_val[13]), 5e-10);
  EXPECT_NEAR(-0.0673880572999235, PetscRealPart(exp_val[14]), 5e-10); 
  EXPECT_NEAR(0.0712550951028285, PetscRealPart(exp_val[15]), 5e-10);
  
  // for (PetscInt i = 0; i < b.nspins; ++i)
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

