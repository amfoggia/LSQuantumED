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

TEST(LatticeTest, NumNeighbours) {
  Environment env{_argc,_argv,12,help};
  
  chain1D l1{env};
  square2D l2{env,4,3};
  square2D l3{env,6,2};
  chain1D l4{env};

  EXPECT_EQ(1,l1.num_neighbours());
  EXPECT_EQ(2,l2.num_neighbours());
  EXPECT_EQ(2,l3.num_neighbours());
  EXPECT_EQ(1,l4.num_neighbours());
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
