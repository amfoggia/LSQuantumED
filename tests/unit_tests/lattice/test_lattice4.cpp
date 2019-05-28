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

TEST(LatticeTest, GetNeighbours1dchain1) {
  Environment env{_argc,_argv,6,help};
  
  chain1D l1{env};
  EXPECT_EQ(12,l1.nn.size());
  
  EXPECT_EQ(0,l1.nn[0]);
  EXPECT_EQ(1,l1.nn[1]);

  EXPECT_EQ(1,l1.nn[2]);
  EXPECT_EQ(2,l1.nn[3]);

  EXPECT_EQ(2,l1.nn[4]);
  EXPECT_EQ(3,l1.nn[5]);

  EXPECT_EQ(3,l1.nn[6]);
  EXPECT_EQ(4,l1.nn[7]);

  EXPECT_EQ(4,l1.nn[8]);
  EXPECT_EQ(5,l1.nn[9]);

  EXPECT_EQ(5,l1.nn[10]);
  EXPECT_EQ(0,l1.nn[11]);

  // Distance vector
  EXPECT_EQ(6,l1.r.size());
  
  EXPECT_EQ(0,l1.r[0][0]);
  EXPECT_EQ(1,l1.r[1][0]);
  EXPECT_EQ(2,l1.r[2][0]);
  EXPECT_EQ(3,l1.r[3][0]);
  EXPECT_EQ(4,l1.r[4][0]);
  EXPECT_EQ(5,l1.r[5][0]);

  EXPECT_EQ(0,l1.r[0][1]);
  EXPECT_EQ(0,l1.r[1][1]);
  EXPECT_EQ(0,l1.r[2][1]);
  EXPECT_EQ(0,l1.r[3][1]);
  EXPECT_EQ(0,l1.r[4][1]);
  EXPECT_EQ(0,l1.r[5][1]);

  // Lattice vector
  EXPECT_EQ(6,l1.q.size());
  
  EXPECT_NEAR(0.0,l1.get_q(0)[0],1e-13);
  EXPECT_NEAR(1.0471975511965976,l1.get_q(1)[0],1e-13);
  EXPECT_NEAR(2.0943951023931953,l1.get_q(2)[0],1e-13);
  EXPECT_NEAR(3.1415926535897931,l1.get_q(3)[0],1e-13);
  EXPECT_NEAR(4.1887902047863905,l1.get_q(4)[0],1e-13);
  EXPECT_NEAR(5.2359877559829879,l1.get_q(5)[0],1e-13);

  EXPECT_NEAR(0,l1.get_q(0)[1],1e-13);
  EXPECT_NEAR(0,l1.get_q(1)[1],1e-13);
  EXPECT_NEAR(0,l1.get_q(2)[1],1e-13);
  EXPECT_NEAR(0,l1.get_q(3)[1],1e-13);
  EXPECT_NEAR(0,l1.get_q(4)[1],1e-13);
  EXPECT_NEAR(0,l1.get_q(5)[1],1e-13);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
