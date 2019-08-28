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

TEST(SolverDiagonal, square2D_J1D1only) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  env.nspins = 16;
  Basis b{env,0};
  square2D lat{env,4,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
  ierr = EPSSetDimensions(solver,2,PETSC_DECIDE,PETSC_DECIDE);

  ierr = Solver::solve(env,solver,nconv); ASSERT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); ASSERT_EQ(0,ierr);

  ASSERT_GT(nconv,1);
  if (nconv > 1) {
    
    ierr = Solver::solution(solver,0,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-11.228483208428859, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver,1,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-10.649884872663442, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
}

TEST(SolverDiagonal, square2D_mag2) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 1.0;
  PetscReal Delta2 = 1.0;
  env.nspins = 16;
  Basis b{env,2};
  square2D lat{env,4,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);
    
  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);

  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
  ierr = EPSSetDimensions(solver,2,PETSC_DECIDE,PETSC_DECIDE);

  ierr = Solver::solve(env,solver,nconv); ASSERT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); ASSERT_EQ(0,ierr);

  ASSERT_GT(nconv,1);
  if (nconv > 1) {
    
    ierr = Solver::solution(solver,0,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-9.8063861588619439, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver,1,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-9.4638734529220887, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
}

TEST(SolverDiagonal, square2D_mag0) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.0;
  PetscReal J2 = 1.0;
  PetscReal Delta2 = 1.0;
  env.nspins = 16;
  Basis b{env,0};
  square2D lat{env,4,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);
  
  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
  ierr = EPSSetDimensions(solver,2,PETSC_DECIDE,PETSC_DECIDE);

  ierr = Solver::solve(env,solver,nconv); ASSERT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); ASSERT_EQ(0,ierr);

  ASSERT_GT(nconv,1);
  if (nconv > 1) {
    
    ierr = Solver::solution(solver,0,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-12.295490216102477, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver,1,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-11.608806942981085, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
}

TEST(SolverDiagonal, square2D_1) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = -0.43;
  PetscReal J2 = -1.0;
  PetscReal Delta2 = 1.0;
  env.nspins = 16;
  Basis b{env,0};
  square2D lat{env,4,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
  ierr = EPSSetDimensions(solver,2,PETSC_DECIDE,PETSC_DECIDE);

  ierr = Solver::solve(env,solver,nconv); ASSERT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); ASSERT_EQ(0,ierr);

  ASSERT_GT(nconv,1);
  if (nconv > 1) {
    
    ierr = Solver::solution(solver,0,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-18.143052224709869, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver,1,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-14.655898142307613, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
}


int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&_argc, _argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}
