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


TEST(LatticeTest, GetType) {
  Environment env{_argc,_argv,12,help};
  
  chain1D l1{env};
  square2D l2{env,3,4};
  square2D l3{env,3,4};
  chain1D l4{env};

  EXPECT_EQ(lattice_type::chain1D, l1.get_type());
  EXPECT_EQ(lattice_type::square2D, l2.get_type());
  EXPECT_EQ(lattice_type::square2D, l3.get_type());
  EXPECT_EQ(lattice_type::chain1D, l4.get_type());
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
