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

// TEST(SolverOffdiag, square2D_lx1_ly6) {
//   PetscReal J1 = 1.0;
//   PetscReal Delta1 = 1.0;
//   PetscReal J2 = 0.0;
//   PetscReal Delta2 = 0.0;
//   env.nspins = 6;
//   Basis b{env,0};
//   square2D lat{env,6,1};
//   Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
//   PetscErrorCode ierr = 0;
//   EPS solver;
//   PetscInt nconv;
//   PetscScalar real_eig;
//   PetscReal error;

//   ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
//   ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

//   EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
//   ierr = EPSSetDimensions(solver, 10, 12, 12); EXPECT_EQ(0,ierr);
 
//   ierr = Solver::solve(env,solver,nconv); EXPECT_EQ(0,ierr);

//   Vec xr;
//   ierr = MatCreateVecs(h.hamilt,NULL,&xr); EXPECT_EQ(0,ierr);

//   ASSERT_GT(nconv,2);
//   if (nconv > 2) {

//     ierr = Solver::solution(solver, 0, real_eig, xr, error); EXPECT_EQ(0,ierr);
//     EXPECT_NEAR(-2.780151527827939, PetscRealPart(real_eig), 1e-13);
//     EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

//     ierr = Solver::solution(solver, 1, real_eig, xr, error); EXPECT_EQ(0,ierr);
//     EXPECT_NEAR(-2.7801515278279414, PetscRealPart(real_eig), 1e-13);
//     EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

//     ierr = Solver::solution(solver, 2, real_eig, xr, error); EXPECT_EQ(0,ierr);
//     EXPECT_NEAR(-2.5223660704270268, PetscRealPart(real_eig), 1e-13);
//     EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

//   }

//   EXPECT_EQ(0, ierr);
//   EXPECT_EQ(0,Solver::SolverClean(solver));
//   EXPECT_EQ(0,VecDestroy(&xr));
// }

// TEST(SolverDelta, square2D_lx6_ly1) {
//   PetscReal J1 = 1.0;
//   PetscReal Delta1 = 1.0;
//   PetscReal J2 = 0.0;
//   PetscReal Delta2 = 0.0;
//   env.nspins = 6;
//   Basis b{env,0};
//   square2D lat{env,6,1};
//   Hamiltonian<square2D> h{env,b,lat,J1,Delta1,J2,Delta2};
//   PetscErrorCode ierr = 0;
//   EPS solver;
//   PetscInt nconv;
//   PetscScalar real_eig;
//   PetscReal error;

//   ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
//   ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

//   EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
//   ierr = EPSSetDimensions(solver, 10, 12, 12); EXPECT_EQ(0,ierr);

//   ierr = Solver::solve(env,solver,nconv); EXPECT_EQ(0,ierr);

//   Vec xr;
//   ierr = MatCreateVecs(h.hamilt,NULL,&xr); EXPECT_EQ(0,ierr);

//   ASSERT_GT(nconv,2);
//   if (nconv > 2) {

//     ierr = Solver::solution(solver, 0, real_eig, xr, error); EXPECT_EQ(0,ierr);
//     EXPECT_NEAR(-2.0072438526146272, PetscRealPart(real_eig), 1e-13);
//     EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

//     ierr = Solver::solution(solver, 1, real_eig, xr, error); EXPECT_EQ(0,ierr);
//     EXPECT_NEAR(-2.0072438526146277, PetscRealPart(real_eig), 1e-13);
//     EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);

//     ierr = Solver::solution(solver, 2, real_eig, xr, error); EXPECT_EQ(0,ierr);
//     EXPECT_NEAR(-1.7195965354971541, PetscRealPart(real_eig), 1e-13);
//     EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
//   }

//   EXPECT_EQ(0, ierr);
//   EXPECT_EQ(0,Solver::SolverClean(solver));
//   EXPECT_EQ(0,VecDestroy(&xr));
// }

TEST(SolverDiagonal, chain1D) {
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
  PetscInt nconv;
  PetscScalar real_eig;
  PetscReal error;
  PetscReal * ref_eigs;

  ierr = PetscCalloc1(9, &ref_eigs); ASSERT_EQ(0,ierr);
  ref_eigs[0] = -2.80277563773199434;
  ref_eigs[1] = -2.11803398874989490;
  ref_eigs[2] = -1.50000000000000067;
  ref_eigs[3] = -1.28077640640441559;
  ref_eigs[4] = -1.28077640640441492;
  ref_eigs[5] = -1.00000000000000044;
  ref_eigs[6] = -0.99999999999999900;
  ref_eigs[7] = -0.50000000000000033;
  ref_eigs[8] = -0.50000000000000011;
  ref_eigs[9] = -0.49999999999999994;
  
  ierr = h.build_diag(env); ASSERT_EQ(0,ierr);
  ierr = h.build_off_diag(env); ASSERT_EQ(0,ierr);

  EXPECT_EQ(0,Solver::SolverInit(solver,h.hamilt,3));
  ierr = EPSSetDimensions(solver, 10, 12, 12); EXPECT_EQ(0,ierr);

  ierr = Solver::solve(env,solver,nconv); EXPECT_EQ(0,ierr);
  
  Vec xr;
  ierr = MatCreateVecs(h.hamilt,NULL,&xr); ASSERT_EQ(0,ierr);

  ASSERT_GT(nconv,1);
  for (int i = 0; i < nconv; ++i) {
    ierr = Solver::solution(solver,i,real_eig,xr,error); ASSERT_EQ(0,ierr);
    EXPECT_NEAR(ref_eigs[i], PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0,Solver::SolverClean(solver));
  EXPECT_EQ(0,VecDestroy(&xr));
  EXPECT_EQ(0,PetscFree(ref_eigs));
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&_argc, _argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}
