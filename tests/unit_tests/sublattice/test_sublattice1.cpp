#include <gtest/gtest.h>
#include "environment.hpp"
#include "sublattice.hpp"

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Sublattice class\n\n";
int _argc;
char ** _argv;

TEST(SublatticeTest, GetSublatsChain1D) {
  Environment env{_argc,_argv,6,help};
  
  chain1D l1{env};
  AF<chain1D> sublat{l1,3};
  EXPECT_EQ(2,sublat.get_size());
  
  EXPECT_EQ(0,sublat.get_sl(0)[0]);
  EXPECT_EQ(2,sublat.get_sl(0)[1]);
  EXPECT_EQ(4,sublat.get_sl(0)[2]);
  
  EXPECT_EQ(1,sublat.get_sl(1)[0]);
  EXPECT_EQ(3,sublat.get_sl(1)[1]);
  EXPECT_EQ(5,sublat.get_sl(1)[2]);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
