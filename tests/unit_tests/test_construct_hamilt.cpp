#include <gtest/gtest.h>
#include "hamiltonian.hpp"
#include <slepcsys.h>

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

TEST(HamiltonianTest, Constructor) {
  Basis b{env,0};
  chain1D lat{env};
  PetscReal j1{1.0}, delta1{0.0}, j2{0.0}, delta2{0.0};

  EXPECT_NO_THROW(Hamiltonian<chain1D> h(env,b,lat,j1,delta1,j2,delta2));
}

TEST(HamiltonianTest, RetrieveAttributes) {
  Basis b{env,0};
  chain1D lat{env};
  PetscReal j1{1.0}, delta1{0.0}, j2{0.0}, delta2{0.0};
  Hamiltonian<chain1D> h{env,b,lat,j1,delta1,j2,delta2};

  EXPECT_EQ(2, h.mpi_size); // In meson.build we run with 2 procs
  if (h.mpi_rank == 0)
    EXPECT_EQ(0, h.mpi_rank);
  else
    EXPECT_EQ(1, h.mpi_rank);

  EXPECT_EQ(20,h.hamilt_size());
  EXPECT_EQ(400,h.hamilt_dim());
  EXPECT_EQ(lattice_type::chain1D,h.hamilt_type());
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}
