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

TEST(LatticeTest, GetNeighbours2Dsquare1) {
  Environment env{_argc,_argv,6,help};
  
  square2D l2{env,6,1};

  // Nearest neighbours
  EXPECT_EQ(18,l2.nn.size());
  
  EXPECT_EQ(0,l2.nn[0]);
  EXPECT_EQ(0,l2.nn[1]);
  EXPECT_EQ(1,l2.nn[2]);

  EXPECT_EQ(1,l2.nn[3]);
  EXPECT_EQ(1,l2.nn[4]);
  EXPECT_EQ(2,l2.nn[5]);

  EXPECT_EQ(2,l2.nn[6]);
  EXPECT_EQ(2,l2.nn[7]);
  EXPECT_EQ(3,l2.nn[8]);

  EXPECT_EQ(3,l2.nn[9]);
  EXPECT_EQ(3,l2.nn[10]);
  EXPECT_EQ(4,l2.nn[11]);

  EXPECT_EQ(4,l2.nn[12]);
  EXPECT_EQ(4,l2.nn[13]);
  EXPECT_EQ(5,l2.nn[14]);

  EXPECT_EQ(5,l2.nn[15]);
  EXPECT_EQ(5,l2.nn[16]);
  EXPECT_EQ(0,l2.nn[17]);

  // Next-nearest neighbours
  EXPECT_EQ(18,l2.nnn.size());
  
  EXPECT_EQ(0,l2.nnn[0]);
  EXPECT_EQ(0,l2.nnn[1]);
  EXPECT_EQ(2,l2.nnn[2]);

  EXPECT_EQ(1,l2.nnn[3]);
  EXPECT_EQ(1,l2.nnn[4]);
  EXPECT_EQ(3,l2.nnn[5]);

  EXPECT_EQ(2,l2.nnn[6]);
  EXPECT_EQ(2,l2.nnn[7]);
  EXPECT_EQ(4,l2.nnn[8]);

  EXPECT_EQ(3,l2.nnn[9]);
  EXPECT_EQ(3,l2.nnn[10]);
  EXPECT_EQ(5,l2.nnn[11]);

  EXPECT_EQ(4,l2.nnn[12]);
  EXPECT_EQ(4,l2.nnn[13]);
  EXPECT_EQ(0,l2.nnn[14]);

  EXPECT_EQ(5,l2.nnn[15]);
  EXPECT_EQ(5,l2.nnn[16]);
  EXPECT_EQ(1,l2.nnn[17]);

  // Distance vector
  // EXPECT_EQ(30,l2.r.size());
  
  // EXPECT_EQ(0,l2.r[0][0]);
  // EXPECT_EQ(0,l2.r[1][0]);
  // EXPECT_EQ(0,l2.r[2][0]);
  // EXPECT_EQ(0,l2.r[3][0]);
  // EXPECT_EQ(0,l2.r[4][0]);
  // EXPECT_EQ(1,l2.r[5][0]);
  // EXPECT_EQ(1,l2.r[6][0]);
  // EXPECT_EQ(1,l2.r[7][0]);
  // EXPECT_EQ(1,l2.r[8][0]);
  // EXPECT_EQ(1,l2.r[9][0]);
  // EXPECT_EQ(2,l2.r[10][0]);
  // EXPECT_EQ(2,l2.r[11][0]);
  // EXPECT_EQ(2,l2.r[12][0]);
  // EXPECT_EQ(2,l2.r[13][0]);
  // EXPECT_EQ(2,l2.r[14][0]);
  // EXPECT_EQ(3,l2.r[15][0]);
  // EXPECT_EQ(3,l2.r[16][0]);
  // EXPECT_EQ(3,l2.r[17][0]);
  // EXPECT_EQ(3,l2.r[18][0]);
  // EXPECT_EQ(3,l2.r[19][0]);
  // EXPECT_EQ(4,l2.r[20][0]);
  // EXPECT_EQ(4,l2.r[21][0]);
  // EXPECT_EQ(4,l2.r[22][0]);
  // EXPECT_EQ(4,l2.r[23][0]);
  // EXPECT_EQ(4,l2.r[24][0]);
  // EXPECT_EQ(5,l2.r[25][0]);
  // EXPECT_EQ(5,l2.r[26][0]);
  // EXPECT_EQ(5,l2.r[27][0]);
  // EXPECT_EQ(5,l2.r[28][0]);
  // EXPECT_EQ(5,l2.r[29][0]);

  // EXPECT_EQ(0,l2.r[0][1]);
  // EXPECT_EQ(1,l2.r[1][1]);
  // EXPECT_EQ(2,l2.r[2][1]);
  // EXPECT_EQ(3,l2.r[3][1]);
  // EXPECT_EQ(4,l2.r[4][1]);
  // EXPECT_EQ(0,l2.r[5][1]);
  // EXPECT_EQ(1,l2.r[6][1]);
  // EXPECT_EQ(2,l2.r[7][1]);
  // EXPECT_EQ(3,l2.r[8][1]);
  // EXPECT_EQ(4,l2.r[9][1]);
  // EXPECT_EQ(0,l2.r[10][1]);
  // EXPECT_EQ(1,l2.r[11][1]);
  // EXPECT_EQ(2,l2.r[12][1]);
  // EXPECT_EQ(3,l2.r[13][1]);
  // EXPECT_EQ(4,l2.r[14][1]);
  // EXPECT_EQ(0,l2.r[15][1]);
  // EXPECT_EQ(1,l2.r[16][1]);
  // EXPECT_EQ(2,l2.r[17][1]);
  // EXPECT_EQ(3,l2.r[18][1]);
  // EXPECT_EQ(4,l2.r[19][1]);
  // EXPECT_EQ(0,l2.r[20][1]);
  // EXPECT_EQ(1,l2.r[21][1]);
  // EXPECT_EQ(2,l2.r[22][1]);
  // EXPECT_EQ(3,l2.r[23][1]);
  // EXPECT_EQ(4,l2.r[24][1]);
  // EXPECT_EQ(0,l2.r[25][1]);
  // EXPECT_EQ(1,l2.r[26][1]);
  // EXPECT_EQ(2,l2.r[27][1]);
  // EXPECT_EQ(3,l2.r[28][1]);
  // EXPECT_EQ(4,l2.r[29][1]);

  // // Lattice vector
  // EXPECT_EQ(30,l2.q.size());
  
  // EXPECT_NEAR(0.0,l2.get_q(0)[0],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(1)[0],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(2)[0],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(3)[0],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(4)[0],1e-13);
  // EXPECT_NEAR(1.0471975511965976,l2.get_q(5)[0],1e-13);
  // EXPECT_NEAR(1.0471975511965976,l2.get_q(6)[0],1e-13);
  // EXPECT_NEAR(1.0471975511965976,l2.get_q(7)[0],1e-13);
  // EXPECT_NEAR(1.0471975511965976,l2.get_q(8)[0],1e-13);
  // EXPECT_NEAR(1.0471975511965976,l2.get_q(9)[0],1e-13);
  // EXPECT_NEAR(2.0943951023931953,l2.get_q(10)[0],1e-13);
  // EXPECT_NEAR(2.0943951023931953,l2.get_q(11)[0],1e-13);
  // EXPECT_NEAR(2.0943951023931953,l2.get_q(12)[0],1e-13);
  // EXPECT_NEAR(2.0943951023931953,l2.get_q(13)[0],1e-13);
  // EXPECT_NEAR(2.0943951023931953,l2.get_q(14)[0],1e-13);
  // EXPECT_NEAR(3.1415926535897931,l2.get_q(15)[0],1e-13);
  // EXPECT_NEAR(3.1415926535897931,l2.get_q(16)[0],1e-13);
  // EXPECT_NEAR(3.1415926535897931,l2.get_q(17)[0],1e-13);
  // EXPECT_NEAR(3.1415926535897931,l2.get_q(18)[0],1e-13);
  // EXPECT_NEAR(3.1415926535897931,l2.get_q(19)[0],1e-13);
  // EXPECT_NEAR(4.1887902047863905,l2.get_q(20)[0],1e-13);
  // EXPECT_NEAR(4.1887902047863905,l2.get_q(21)[0],1e-13);
  // EXPECT_NEAR(4.1887902047863905,l2.get_q(22)[0],1e-13);
  // EXPECT_NEAR(4.1887902047863905,l2.get_q(23)[0],1e-13);
  // EXPECT_NEAR(4.1887902047863905,l2.get_q(24)[0],1e-13);
  // EXPECT_NEAR(5.2359877559829879,l2.get_q(25)[0],1e-13);
  // EXPECT_NEAR(5.2359877559829879,l2.get_q(26)[0],1e-13);
  // EXPECT_NEAR(5.2359877559829879,l2.get_q(27)[0],1e-13);
  // EXPECT_NEAR(5.2359877559829879,l2.get_q(28)[0],1e-13);
  // EXPECT_NEAR(5.2359877559829879,l2.get_q(29)[0],1e-13);

  // EXPECT_NEAR(0.0,l2.get_q(0)[1],1e-13);
  // EXPECT_NEAR(1.2566370614359172,l2.get_q(1)[1],1e-13);
  // EXPECT_NEAR(2.5132741228718345,l2.get_q(2)[1],1e-13);
  // EXPECT_NEAR(3.7699111843077517,l2.get_q(3)[1],1e-13);
  // EXPECT_NEAR(5.026548245743669,l2.get_q(4)[1],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(5)[1],1e-13);
  // EXPECT_NEAR(1.2566370614359172,l2.get_q(6)[1],1e-13);
  // EXPECT_NEAR(2.5132741228718345,l2.get_q(7)[1],1e-13);
  // EXPECT_NEAR(3.7699111843077517,l2.get_q(8)[1],1e-13);
  // EXPECT_NEAR(5.026548245743669,l2.get_q(9)[1],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(10)[1],1e-13);
  // EXPECT_NEAR(1.2566370614359172,l2.get_q(11)[1],1e-13);
  // EXPECT_NEAR(2.5132741228718345,l2.get_q(12)[1],1e-13);
  // EXPECT_NEAR(3.7699111843077517,l2.get_q(13)[1],1e-13);
  // EXPECT_NEAR(5.026548245743669,l2.get_q(14)[1],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(15)[1],1e-13);
  // EXPECT_NEAR(1.2566370614359172,l2.get_q(16)[1],1e-13);
  // EXPECT_NEAR(2.5132741228718345,l2.get_q(17)[1],1e-13);
  // EXPECT_NEAR(3.7699111843077517,l2.get_q(18)[1],1e-13);
  // EXPECT_NEAR(5.026548245743669,l2.get_q(19)[1],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(20)[1],1e-13);
  // EXPECT_NEAR(1.2566370614359172,l2.get_q(21)[1],1e-13);
  // EXPECT_NEAR(2.5132741228718345,l2.get_q(22)[1],1e-13);
  // EXPECT_NEAR(3.7699111843077517,l2.get_q(23)[1],1e-13);
  // EXPECT_NEAR(5.026548245743669,l2.get_q(24)[1],1e-13);
  // EXPECT_NEAR(0.0,l2.get_q(25)[1],1e-13);
  // EXPECT_NEAR(1.2566370614359172,l2.get_q(26)[1],1e-13);
  // EXPECT_NEAR(2.5132741228718345,l2.get_q(27)[1],1e-13);
  // EXPECT_NEAR(3.7699111843077517,l2.get_q(28)[1],1e-13);
  // EXPECT_NEAR(5.026548245743669,l2.get_q(29)[1],1e-13);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
