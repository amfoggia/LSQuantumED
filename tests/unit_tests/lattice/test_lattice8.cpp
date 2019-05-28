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

TEST(LatticeTest, GetNNN1dchain1) {
  Environment env{_argc,_argv,8,help};
  
  chain1D l1{env};
  ASSERT_EQ(16,l1.nnn.size());

  EXPECT_EQ(0,l1.nnn[0]);
  EXPECT_EQ(2,l1.nnn[1]);

  EXPECT_EQ(1,l1.nnn[2]);
  EXPECT_EQ(3,l1.nnn[3]);

  EXPECT_EQ(2,l1.nnn[4]);
  EXPECT_EQ(4,l1.nnn[5]);

  EXPECT_EQ(3,l1.nnn[6]);
  EXPECT_EQ(5,l1.nnn[7]);

  EXPECT_EQ(4,l1.nnn[8]);
  EXPECT_EQ(6,l1.nnn[9]);

  EXPECT_EQ(5,l1.nnn[10]);
  EXPECT_EQ(7,l1.nnn[11]);

  EXPECT_EQ(6,l1.nnn[12]);
  EXPECT_EQ(0,l1.nnn[13]);

  EXPECT_EQ(7,l1.nnn[14]);
  EXPECT_EQ(1,l1.nnn[15]);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
