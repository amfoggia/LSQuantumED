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
  PetscScalar real_eig;
  PetscReal error;

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
    EXPECT_NEAR(-2.8027756377319952, PetscRealPart(real_eig), 1e-13);
    EXPECT_PRED_FORMAT2(::testing::FloatLE, (double)error, 1e-8);
  }

  std::vector<PetscInt> tmp;
  for (PetscInt i = 0; i < env.nspins; ++i)
    tmp.push_back(i);
  
  PetscScalar dynstructfactor;
  Phys::DSF_data<chain1D> dsf_data{};
  dsf_data.b0 = &b;
  dsf_data.b1 = &b;
  dsf_data.lat = &lat;
  dsf_data.nQ = 6;
  dsf_data.qi = &tmp;
  dsf_data.nev = 10;
  dsf_data.ncv = PETSC_DECIDE;
  dsf_data.mpd = PETSC_DECIDE;
  dsf_data.h0 = &h;
  dsf_data.h1 = &h;
  dsf_data.disorder = PETSC_FALSE;
  dsf_data.hi = NULL;
  dsf_data.path = "./";
  
  dynstructfactor = Phys::DynStructFactor<chain1D,Sqz>(env, dsf_data, real_eig, xr, 0);
  
  EXPECT_EQ(0, ierr);
  EXPECT_EQ(0, Solver::SolverClean(solver));
  EXPECT_EQ(0, VecDestroy(&xr));
  
}


int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

