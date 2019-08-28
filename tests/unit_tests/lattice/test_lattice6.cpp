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
  EXPECT_EQ(30,l2.nnXsite.size());
  
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

  for (PetscInt i = 0; i < env.nspins; ++i)
    EXPECT_EQ(2,l2.num_nnXsite(i));

  // Next-nearest neighbours
  EXPECT_EQ(90,l2.nnn.size());
  EXPECT_EQ(30,l2.nnnXsite.size());
  
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

  for (PetscInt i = 0; i < env.nspins; ++i)
    EXPECT_EQ(2,l2.num_nnnXsite(i));
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
