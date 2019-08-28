#include "hamiltonian.hpp"
#include "solver.hpp"
#include "phys.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Correlation-SzSz class\n\n";
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


TEST(DSF, chain1D) {
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
  PetscScalar real_eig, im_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

  MatView(h.hamilt, PETSC_VIEWER_STDOUT_WORLD);
  
  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,1));
 
  ierr = Solver::solve(env,solver,nconv); EXPECT_EQ(0,ierr);

  Vec xr, xi;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); EXPECT_EQ(0,ierr);
  ierr = MatCreateVecs(h.hamilt,NULL,&xi); EXPECT_EQ(0,ierr);

  ASSERT_GE(nconv,1);
  if (nconv >= 1) {
    //    ierr = Solver::solution(solver, 0, real_eig, xr, error); EXPECT_EQ(0,ierr);
    ierr = EPSGetEigenpair(solver, 0, &real_eig, &im_eig, xr, xi);
    EXPECT_NEAR(-2.8027756377319952, PetscRealPart(real_eig), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(real_eig), 1e-13);
    PetscPrintf(PETSC_COMM_WORLD, "real_eig: %f, im_eig: %f\n", real_eig, im_eig);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
    //ierr = VecView(xr, PETSC_VIEWER_STDOUT_WORLD);
    //ierr = VecView(xi, PETSC_VIEWER_STDOUT_WORLD);
    //VecDot(xr, xr, &vecdot);
    //PetscPrintf(PETSC_COMM_WORLD, "vec_dot: %f\n", vecdot);
  }
  
  Phys::DSF_data<chain1D> dsf_data{};
  dsf_data.b0 = &b;
  dsf_data.b1 = &b;
  dsf_data.lat = &lat;
  dsf_data.nev = 10;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
  dsf_data.h0 = &h;
  dsf_data.h1 = &h;
  dsf_data.disorder = PETSC_FALSE;
  dsf_data.hi = NULL;
  dsf_data.path = "./dsf1";
  
  ierr = Phys::DynStructFactor<chain1D,Sqz>(env, dsf_data, real_eig, xr, 0);
  
  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0, Solver::SolverClean(solver));
  EXPECT_EQ(0, VecDestroy(&xr));
  EXPECT_EQ(0, VecDestroy(&xi));
  
}


int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

