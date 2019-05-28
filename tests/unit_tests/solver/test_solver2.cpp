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

Environment env{_argc,_argv,12,help};

class HamiltonianTestEnv : public ::testing::Environment {
protected:
  virtual void TearDown() {
    ::testing::internal::CaptureStdout();
    EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};

TEST(SolverOffdiag, square2Ddiag0) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 0.0;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  Basis b{env,4};
  square2D lat{env,3,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

#ifdef DISORDER
  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
#endif
  
  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
  ierr = EPSSetDimensions(solver, 10, 12, 12); EXPECT_EQ(0,ierr);
 
  ierr = Solver::solve_lanczos(env,solver,nconv); EXPECT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); EXPECT_EQ(0,ierr);

  ASSERT_GT(nconv,2);
  if (nconv > 2) {

    ierr = Solver::solution(solver, 0, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-2.780151527827939, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver, 1, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-2.7801515278279414, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver, 2, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-2.5223660704270268, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
}

TEST(SolverDelta, square2Dfull) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 0.3456;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  Basis b{env,4};
  square2D lat{env,3,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

#ifdef DISORDER
  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
#endif

  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
  ierr = EPSSetDimensions(solver, 10, 12, 12); EXPECT_EQ(0,ierr);

  ierr = Solver::solve_lanczos(env,solver,nconv); EXPECT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); EXPECT_EQ(0,ierr);

  ASSERT_GT(nconv,2);
  if (nconv > 2) {

    ierr = Solver::solution(solver, 0, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-2.0072438526146272, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver, 1, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-2.0072438526146277, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver, 2, real_eig, xr, error); EXPECT_EQ(0,ierr);
    EXPECT_NEAR(-1.7195965354971541, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
}

#ifdef BUILD_SEPARATE
TEST(SolverDiagonal, square2DNoOffdiag) {
  PetscReal J1 = 1.0;
  PetscReal Delta1 = 1.4;
  PetscReal J2 = 0.0;
  PetscReal Delta2 = 0.0;
  Basis b{env,0};
  square2D lat{env,3,4};
  Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
  PetscErrorCode ierr = 0;
  EPS solver;
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;

  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  
  ierr = MatAssemblyBegin(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);
  ierr = MatAssemblyEnd(h.hamilt,MAT_FINAL_ASSEMBLY); ASSERT_EQ(0,ierr);

  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,2));

  ierr = Solver::solve_lanczos(env,solver,nconv); ASSERT_EQ(0,ierr);

  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); ASSERT_EQ(0,ierr);

  ASSERT_GT(nconv,1);
  if (nconv > 2) {
    
    ierr = Solver::solution(solver,0,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-5.6, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

    ierr = Solver::solution(solver,1,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(-4.2, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
}
#endif

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&_argc, _argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}
