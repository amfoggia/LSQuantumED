#include "environment.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Environment class\n\n";
int _argc;
char ** _argv;

TEST(EnvironmentTest, ConstrDestr84) {
  Environment testEnv(_argc, _argv, 8, help);

  EXPECT_EQ(2, testEnv.mpi_size); // In meson.build we run with 2 procs
  if (testEnv.mpi_rank == 0)
    EXPECT_EQ(0, testEnv.mpi_rank);
  else
    EXPECT_EQ(1, testEnv.mpi_rank);

  EXPECT_EQ(8, testEnv.nspins);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
