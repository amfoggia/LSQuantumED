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
  Environment env{_argc,_argv,30,help};

  square2D l2{env,6,5};

  // Nearest neighbours
  EXPECT_EQ(90,l2.nn.size());
  
  EXPECT_EQ(0,l2.nn[0]);
  EXPECT_EQ(1,l2.nn[1]);
  EXPECT_EQ(5,l2.nn[2]);

  EXPECT_EQ(1,l2.nn[3]);
  EXPECT_EQ(2,l2.nn[4]);
  EXPECT_EQ(6,l2.nn[5]);

  EXPECT_EQ(2,l2.nn[6]);
  EXPECT_EQ(3,l2.nn[7]);
  EXPECT_EQ(7,l2.nn[8]);

  EXPECT_EQ(3,l2.nn[9]);
  EXPECT_EQ(4,l2.nn[10]);
  EXPECT_EQ(8,l2.nn[11]);

  EXPECT_EQ(4,l2.nn[12]);
  EXPECT_EQ(0,l2.nn[13]);
  EXPECT_EQ(9,l2.nn[14]);

  EXPECT_EQ(5,l2.nn[15]);
  EXPECT_EQ(6,l2.nn[16]);
  EXPECT_EQ(10,l2.nn[17]);

  EXPECT_EQ(6,l2.nn[18]);
  EXPECT_EQ(7,l2.nn[19]);
  EXPECT_EQ(11,l2.nn[20]);

  EXPECT_EQ(7,l2.nn[21]);
  EXPECT_EQ(8,l2.nn[22]);
  EXPECT_EQ(12,l2.nn[23]);

  EXPECT_EQ(8,l2.nn[24]);
  EXPECT_EQ(9,l2.nn[25]);
  EXPECT_EQ(13,l2.nn[26]);

  EXPECT_EQ(9,l2.nn[27]);
  EXPECT_EQ(5,l2.nn[28]);
  EXPECT_EQ(14,l2.nn[29]);

  EXPECT_EQ(10,l2.nn[30]);
  EXPECT_EQ(11,l2.nn[31]);
  EXPECT_EQ(15,l2.nn[32]);

  EXPECT_EQ(11,l2.nn[33]);
  EXPECT_EQ(12,l2.nn[34]);
  EXPECT_EQ(16,l2.nn[35]);

  EXPECT_EQ(12,l2.nn[36]);
  EXPECT_EQ(13,l2.nn[37]);
  EXPECT_EQ(17,l2.nn[38]);

  EXPECT_EQ(13,l2.nn[39]);
  EXPECT_EQ(14,l2.nn[40]);
  EXPECT_EQ(18,l2.nn[41]);

  EXPECT_EQ(14,l2.nn[42]);
  EXPECT_EQ(10,l2.nn[43]);
  EXPECT_EQ(19,l2.nn[44]);

  EXPECT_EQ(15,l2.nn[45]);
  EXPECT_EQ(16,l2.nn[46]);
  EXPECT_EQ(20,l2.nn[47]);

  EXPECT_EQ(16,l2.nn[48]);
  EXPECT_EQ(17,l2.nn[49]);
  EXPECT_EQ(21,l2.nn[50]);

  EXPECT_EQ(17,l2.nn[51]);
  EXPECT_EQ(18,l2.nn[52]);
  EXPECT_EQ(22,l2.nn[53]);

  EXPECT_EQ(18,l2.nn[54]);
  EXPECT_EQ(19,l2.nn[55]);
  EXPECT_EQ(23,l2.nn[56]);

  EXPECT_EQ(19,l2.nn[57]);
  EXPECT_EQ(15,l2.nn[58]);
  EXPECT_EQ(24,l2.nn[59]);

  EXPECT_EQ(20,l2.nn[60]);
  EXPECT_EQ(21,l2.nn[61]);
  EXPECT_EQ(25,l2.nn[62]);

  EXPECT_EQ(21,l2.nn[63]);
  EXPECT_EQ(22,l2.nn[64]);
  EXPECT_EQ(26,l2.nn[65]);

  EXPECT_EQ(22,l2.nn[66]);
  EXPECT_EQ(23,l2.nn[67]);
  EXPECT_EQ(27,l2.nn[68]);

  EXPECT_EQ(23,l2.nn[69]);
  EXPECT_EQ(24,l2.nn[70]);
  EXPECT_EQ(28,l2.nn[71]);

  EXPECT_EQ(24,l2.nn[72]);
  EXPECT_EQ(20,l2.nn[73]);
  EXPECT_EQ(29,l2.nn[74]);

  EXPECT_EQ(25,l2.nn[75]);
  EXPECT_EQ(26,l2.nn[76]);
  EXPECT_EQ(0,l2.nn[77]);

  EXPECT_EQ(26,l2.nn[78]);
  EXPECT_EQ(27,l2.nn[79]);
  EXPECT_EQ(1,l2.nn[80]);

  EXPECT_EQ(27,l2.nn[81]);
  EXPECT_EQ(28,l2.nn[82]);
  EXPECT_EQ(2,l2.nn[83]);

  EXPECT_EQ(28,l2.nn[84]);
  EXPECT_EQ(29,l2.nn[85]);
  EXPECT_EQ(3,l2.nn[86]);

  EXPECT_EQ(29,l2.nn[87]);
  EXPECT_EQ(25,l2.nn[88]);
  EXPECT_EQ(4,l2.nn[89]);

  // Next-nearest neighbours
  EXPECT_EQ(90,l2.nnn.size());
  
  EXPECT_EQ(0,l2.nnn[0]);
  EXPECT_EQ(6,l2.nnn[1]);
  EXPECT_EQ(9,l2.nnn[2]);

  EXPECT_EQ(1,l2.nnn[3]);
  EXPECT_EQ(7,l2.nnn[4]);
  EXPECT_EQ(5,l2.nnn[5]);

  EXPECT_EQ(2,l2.nnn[6]);
  EXPECT_EQ(8,l2.nnn[7]);
  EXPECT_EQ(6,l2.nnn[8]);

  EXPECT_EQ(3,l2.nnn[9]);
  EXPECT_EQ(9,l2.nnn[10]);
  EXPECT_EQ(7,l2.nnn[11]);

  EXPECT_EQ(4,l2.nnn[12]);
  EXPECT_EQ(5,l2.nnn[13]);
  EXPECT_EQ(8,l2.nnn[14]);

  EXPECT_EQ(5,l2.nnn[15]);
  EXPECT_EQ(11,l2.nnn[16]);
  EXPECT_EQ(14,l2.nnn[17]);

  EXPECT_EQ(6,l2.nnn[18]);
  EXPECT_EQ(12,l2.nnn[19]);
  EXPECT_EQ(10,l2.nnn[20]);

  EXPECT_EQ(7,l2.nnn[21]);
  EXPECT_EQ(13,l2.nnn[22]);
  EXPECT_EQ(11,l2.nnn[23]);

  EXPECT_EQ(8,l2.nnn[24]);
  EXPECT_EQ(14,l2.nnn[25]);
  EXPECT_EQ(12,l2.nnn[26]);

  EXPECT_EQ(9,l2.nnn[27]);
  EXPECT_EQ(10,l2.nnn[28]);
  EXPECT_EQ(13,l2.nnn[29]);

  EXPECT_EQ(10,l2.nnn[30]);
  EXPECT_EQ(16,l2.nnn[31]);
  EXPECT_EQ(19,l2.nnn[32]);

  EXPECT_EQ(11,l2.nnn[33]);
  EXPECT_EQ(17,l2.nnn[34]);
  EXPECT_EQ(15,l2.nnn[35]);

  EXPECT_EQ(12,l2.nnn[36]);
  EXPECT_EQ(18,l2.nnn[37]);
  EXPECT_EQ(16,l2.nnn[38]);

  EXPECT_EQ(13,l2.nnn[39]);
  EXPECT_EQ(19,l2.nnn[40]);
  EXPECT_EQ(17,l2.nnn[41]);

  EXPECT_EQ(14,l2.nnn[42]);
  EXPECT_EQ(15,l2.nnn[43]);
  EXPECT_EQ(18,l2.nnn[44]);

  EXPECT_EQ(15,l2.nnn[45]);
  EXPECT_EQ(21,l2.nnn[46]);
  EXPECT_EQ(24,l2.nnn[47]);

  EXPECT_EQ(16,l2.nnn[48]);
  EXPECT_EQ(22,l2.nnn[49]);
  EXPECT_EQ(20,l2.nnn[50]);

  EXPECT_EQ(17,l2.nnn[51]);
  EXPECT_EQ(23,l2.nnn[52]);
  EXPECT_EQ(21,l2.nnn[53]);

  EXPECT_EQ(18,l2.nnn[54]);
  EXPECT_EQ(24,l2.nnn[55]);
  EXPECT_EQ(22,l2.nnn[56]);

  EXPECT_EQ(19,l2.nnn[57]);
  EXPECT_EQ(20,l2.nnn[58]);
  EXPECT_EQ(23,l2.nnn[59]);

  EXPECT_EQ(20,l2.nnn[60]);
  EXPECT_EQ(26,l2.nnn[61]);
  EXPECT_EQ(29,l2.nnn[62]);

  EXPECT_EQ(21,l2.nnn[63]);
  EXPECT_EQ(27,l2.nnn[64]);
  EXPECT_EQ(25,l2.nnn[65]);

  EXPECT_EQ(22,l2.nnn[66]);
  EXPECT_EQ(28,l2.nnn[67]);
  EXPECT_EQ(26,l2.nnn[68]);

  EXPECT_EQ(23,l2.nnn[69]);
  EXPECT_EQ(29,l2.nnn[70]);
  EXPECT_EQ(27,l2.nnn[71]);

  EXPECT_EQ(24,l2.nnn[72]);
  EXPECT_EQ(25,l2.nnn[73]);
  EXPECT_EQ(28,l2.nnn[74]);

  EXPECT_EQ(25,l2.nnn[75]);
  EXPECT_EQ(1,l2.nnn[76]);
  EXPECT_EQ(4,l2.nnn[77]);

  EXPECT_EQ(26,l2.nnn[78]);
  EXPECT_EQ(2,l2.nnn[79]);
  EXPECT_EQ(0,l2.nnn[80]);

  EXPECT_EQ(27,l2.nnn[81]);
  EXPECT_EQ(3,l2.nnn[82]);
  EXPECT_EQ(1,l2.nnn[83]);

  EXPECT_EQ(28,l2.nnn[84]);
  EXPECT_EQ(4,l2.nnn[85]);
  EXPECT_EQ(2,l2.nnn[86]);

  EXPECT_EQ(29,l2.nnn[87]);
  EXPECT_EQ(0,l2.nnn[88]);
  EXPECT_EQ(3,l2.nnn[89]);

  // Distance vector
  EXPECT_EQ(30,l2.r.size());
  
  EXPECT_EQ(0,l2.r[0][0]);
  EXPECT_EQ(0,l2.r[1][0]);
  EXPECT_EQ(0,l2.r[2][0]);
  EXPECT_EQ(0,l2.r[3][0]);
  EXPECT_EQ(0,l2.r[4][0]);
  EXPECT_EQ(1,l2.r[5][0]);
  EXPECT_EQ(1,l2.r[6][0]);
  EXPECT_EQ(1,l2.r[7][0]);
  EXPECT_EQ(1,l2.r[8][0]);
  EXPECT_EQ(1,l2.r[9][0]);
  EXPECT_EQ(2,l2.r[10][0]);
  EXPECT_EQ(2,l2.r[11][0]);
  EXPECT_EQ(2,l2.r[12][0]);
  EXPECT_EQ(2,l2.r[13][0]);
  EXPECT_EQ(2,l2.r[14][0]);
  EXPECT_EQ(3,l2.r[15][0]);
  EXPECT_EQ(3,l2.r[16][0]);
  EXPECT_EQ(3,l2.r[17][0]);
  EXPECT_EQ(3,l2.r[18][0]);
  EXPECT_EQ(3,l2.r[19][0]);
  EXPECT_EQ(4,l2.r[20][0]);
  EXPECT_EQ(4,l2.r[21][0]);
  EXPECT_EQ(4,l2.r[22][0]);
  EXPECT_EQ(4,l2.r[23][0]);
  EXPECT_EQ(4,l2.r[24][0]);
  EXPECT_EQ(5,l2.r[25][0]);
  EXPECT_EQ(5,l2.r[26][0]);
  EXPECT_EQ(5,l2.r[27][0]);
  EXPECT_EQ(5,l2.r[28][0]);
  EXPECT_EQ(5,l2.r[29][0]);

  EXPECT_EQ(0,l2.r[0][1]);
  EXPECT_EQ(1,l2.r[1][1]);
  EXPECT_EQ(2,l2.r[2][1]);
  EXPECT_EQ(3,l2.r[3][1]);
  EXPECT_EQ(4,l2.r[4][1]);
  EXPECT_EQ(0,l2.r[5][1]);
  EXPECT_EQ(1,l2.r[6][1]);
  EXPECT_EQ(2,l2.r[7][1]);
  EXPECT_EQ(3,l2.r[8][1]);
  EXPECT_EQ(4,l2.r[9][1]);
  EXPECT_EQ(0,l2.r[10][1]);
  EXPECT_EQ(1,l2.r[11][1]);
  EXPECT_EQ(2,l2.r[12][1]);
  EXPECT_EQ(3,l2.r[13][1]);
  EXPECT_EQ(4,l2.r[14][1]);
  EXPECT_EQ(0,l2.r[15][1]);
  EXPECT_EQ(1,l2.r[16][1]);
  EXPECT_EQ(2,l2.r[17][1]);
  EXPECT_EQ(3,l2.r[18][1]);
  EXPECT_EQ(4,l2.r[19][1]);
  EXPECT_EQ(0,l2.r[20][1]);
  EXPECT_EQ(1,l2.r[21][1]);
  EXPECT_EQ(2,l2.r[22][1]);
  EXPECT_EQ(3,l2.r[23][1]);
  EXPECT_EQ(4,l2.r[24][1]);
  EXPECT_EQ(0,l2.r[25][1]);
  EXPECT_EQ(1,l2.r[26][1]);
  EXPECT_EQ(2,l2.r[27][1]);
  EXPECT_EQ(3,l2.r[28][1]);
  EXPECT_EQ(4,l2.r[29][1]);

  // Lattice vector
  EXPECT_EQ(30,l2.q.size());
  
  EXPECT_NEAR(0.0,l2.get_q(0)[0],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(1)[0],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(2)[0],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(3)[0],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(4)[0],1e-13);
  EXPECT_NEAR(1.0471975511965976,l2.get_q(5)[0],1e-13);
  EXPECT_NEAR(1.0471975511965976,l2.get_q(6)[0],1e-13);
  EXPECT_NEAR(1.0471975511965976,l2.get_q(7)[0],1e-13);
  EXPECT_NEAR(1.0471975511965976,l2.get_q(8)[0],1e-13);
  EXPECT_NEAR(1.0471975511965976,l2.get_q(9)[0],1e-13);
  EXPECT_NEAR(2.0943951023931953,l2.get_q(10)[0],1e-13);
  EXPECT_NEAR(2.0943951023931953,l2.get_q(11)[0],1e-13);
  EXPECT_NEAR(2.0943951023931953,l2.get_q(12)[0],1e-13);
  EXPECT_NEAR(2.0943951023931953,l2.get_q(13)[0],1e-13);
  EXPECT_NEAR(2.0943951023931953,l2.get_q(14)[0],1e-13);
  EXPECT_NEAR(3.1415926535897931,l2.get_q(15)[0],1e-13);
  EXPECT_NEAR(3.1415926535897931,l2.get_q(16)[0],1e-13);
  EXPECT_NEAR(3.1415926535897931,l2.get_q(17)[0],1e-13);
  EXPECT_NEAR(3.1415926535897931,l2.get_q(18)[0],1e-13);
  EXPECT_NEAR(3.1415926535897931,l2.get_q(19)[0],1e-13);
  EXPECT_NEAR(4.1887902047863905,l2.get_q(20)[0],1e-13);
  EXPECT_NEAR(4.1887902047863905,l2.get_q(21)[0],1e-13);
  EXPECT_NEAR(4.1887902047863905,l2.get_q(22)[0],1e-13);
  EXPECT_NEAR(4.1887902047863905,l2.get_q(23)[0],1e-13);
  EXPECT_NEAR(4.1887902047863905,l2.get_q(24)[0],1e-13);
  EXPECT_NEAR(5.2359877559829879,l2.get_q(25)[0],1e-13);
  EXPECT_NEAR(5.2359877559829879,l2.get_q(26)[0],1e-13);
  EXPECT_NEAR(5.2359877559829879,l2.get_q(27)[0],1e-13);
  EXPECT_NEAR(5.2359877559829879,l2.get_q(28)[0],1e-13);
  EXPECT_NEAR(5.2359877559829879,l2.get_q(29)[0],1e-13);

  EXPECT_NEAR(0.0,l2.get_q(0)[1],1e-13);
  EXPECT_NEAR(1.2566370614359172,l2.get_q(1)[1],1e-13);
  EXPECT_NEAR(2.5132741228718345,l2.get_q(2)[1],1e-13);
  EXPECT_NEAR(3.7699111843077517,l2.get_q(3)[1],1e-13);
  EXPECT_NEAR(5.026548245743669,l2.get_q(4)[1],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(5)[1],1e-13);
  EXPECT_NEAR(1.2566370614359172,l2.get_q(6)[1],1e-13);
  EXPECT_NEAR(2.5132741228718345,l2.get_q(7)[1],1e-13);
  EXPECT_NEAR(3.7699111843077517,l2.get_q(8)[1],1e-13);
  EXPECT_NEAR(5.026548245743669,l2.get_q(9)[1],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(10)[1],1e-13);
  EXPECT_NEAR(1.2566370614359172,l2.get_q(11)[1],1e-13);
  EXPECT_NEAR(2.5132741228718345,l2.get_q(12)[1],1e-13);
  EXPECT_NEAR(3.7699111843077517,l2.get_q(13)[1],1e-13);
  EXPECT_NEAR(5.026548245743669,l2.get_q(14)[1],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(15)[1],1e-13);
  EXPECT_NEAR(1.2566370614359172,l2.get_q(16)[1],1e-13);
  EXPECT_NEAR(2.5132741228718345,l2.get_q(17)[1],1e-13);
  EXPECT_NEAR(3.7699111843077517,l2.get_q(18)[1],1e-13);
  EXPECT_NEAR(5.026548245743669,l2.get_q(19)[1],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(20)[1],1e-13);
  EXPECT_NEAR(1.2566370614359172,l2.get_q(21)[1],1e-13);
  EXPECT_NEAR(2.5132741228718345,l2.get_q(22)[1],1e-13);
  EXPECT_NEAR(3.7699111843077517,l2.get_q(23)[1],1e-13);
  EXPECT_NEAR(5.026548245743669,l2.get_q(24)[1],1e-13);
  EXPECT_NEAR(0.0,l2.get_q(25)[1],1e-13);
  EXPECT_NEAR(1.2566370614359172,l2.get_q(26)[1],1e-13);
  EXPECT_NEAR(2.5132741228718345,l2.get_q(27)[1],1e-13);
  EXPECT_NEAR(3.7699111843077517,l2.get_q(28)[1],1e-13);
  EXPECT_NEAR(5.026548245743669,l2.get_q(29)[1],1e-13);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
