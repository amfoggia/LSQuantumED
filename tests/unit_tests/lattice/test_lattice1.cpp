#include <gtest/gtest.h>
#include "lattice.hpp"

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Lattice class\n\n";
int _argc;
char ** _argv;

TEST(LatticeTest, Constructor) {
  Environment env{_argc,_argv,12,help};
  
  EXPECT_NO_THROW(chain1D l2(env));

  // 2d square constructor
  EXPECT_NO_THROW(square2D l2(env,3,4));
  // This one throws an exception because nspins=12 != 4x4
  ASSERT_THROW(square2D l2(env,4,4), std::invalid_argument);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

