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

TEST(LatticeTest, GetNeighbours1dchain2) {
  Environment env{_argc,_argv,4,help};
  
  chain1D l1{env};

  // Nearest neighbours
  EXPECT_EQ(8,l1.nn.size());
  
  EXPECT_EQ(0,l1.nn[0]);
  EXPECT_EQ(1,l1.nn[1]);

  EXPECT_EQ(1,l1.nn[2]);
  EXPECT_EQ(2,l1.nn[3]);

  EXPECT_EQ(2,l1.nn[4]);
  EXPECT_EQ(3,l1.nn[5]);

  EXPECT_EQ(3,l1.nn[6]);
  EXPECT_EQ(0,l1.nn[7]);

  // Next-nearest neighbours
  ASSERT_EQ(8,l1.nnn.size());

  EXPECT_EQ(0,l1.nnn[0]);
  EXPECT_EQ(2,l1.nnn[1]);

  EXPECT_EQ(1,l1.nnn[2]);
  EXPECT_EQ(3,l1.nnn[3]);

  EXPECT_EQ(2,l1.nnn[4]);
  EXPECT_EQ(0,l1.nnn[5]);
  
  EXPECT_EQ(3,l1.nnn[6]);
  EXPECT_EQ(1,l1.nnn[7]);

  // Distance vector
  ASSERT_EQ(4,l1.r.size());

  EXPECT_EQ(0,l1.r[0][0]);
  EXPECT_EQ(1,l1.r[1][0]);
  EXPECT_EQ(2,l1.r[2][0]);
  EXPECT_EQ(3,l1.r[3][0]);

  EXPECT_EQ(0,l1.r[0][1]);
  EXPECT_EQ(0,l1.r[1][1]);
  EXPECT_EQ(0,l1.r[2][1]);
  EXPECT_EQ(0,l1.r[3][1]);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
