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

TEST(LatticeTest, GetNeighbours2Dsquare2) {
  Environment env{_argc,_argv,12,help};

  square2D  l2{env,3,4};

  // Nearest neighbours
  EXPECT_EQ(36,l2.nn.size());
  EXPECT_EQ(12,l2.nnXsite.size());
  
  EXPECT_EQ(0,l2.nn[0]);
  EXPECT_EQ(1,l2.nn[1]);
  EXPECT_EQ(4,l2.nn[2]);

  EXPECT_EQ(1,l2.nn[3]);
  EXPECT_EQ(2,l2.nn[4]);
  EXPECT_EQ(5,l2.nn[5]);

  EXPECT_EQ(2,l2.nn[6]);
  EXPECT_EQ(3,l2.nn[7]);
  EXPECT_EQ(6,l2.nn[8]);

  EXPECT_EQ(3,l2.nn[9]);
  EXPECT_EQ(0,l2.nn[10]);
  EXPECT_EQ(7,l2.nn[11]);

  EXPECT_EQ(4,l2.nn[12]);
  EXPECT_EQ(5,l2.nn[13]);
  EXPECT_EQ(8,l2.nn[14]);

  EXPECT_EQ(5,l2.nn[15]);
  EXPECT_EQ(6,l2.nn[16]);
  EXPECT_EQ(9,l2.nn[17]);

  EXPECT_EQ(6,l2.nn[18]);
  EXPECT_EQ(7,l2.nn[19]);
  EXPECT_EQ(10,l2.nn[20]);

  EXPECT_EQ(7,l2.nn[21]);
  EXPECT_EQ(4,l2.nn[22]);
  EXPECT_EQ(11,l2.nn[23]);

  EXPECT_EQ(8,l2.nn[24]);
  EXPECT_EQ(9,l2.nn[25]);
  EXPECT_EQ(0,l2.nn[26]);
  
  EXPECT_EQ(9,l2.nn[27]);
  EXPECT_EQ(10,l2.nn[28]);
  EXPECT_EQ(1,l2.nn[29]);

  EXPECT_EQ(10,l2.nn[30]);
  EXPECT_EQ(11,l2.nn[31]);
  EXPECT_EQ(2,l2.nn[32]);

  EXPECT_EQ(11,l2.nn[33]);
  EXPECT_EQ(8,l2.nn[34]);
  EXPECT_EQ(3,l2.nn[35]);

  for (PetscInt i = 0; i < env.nspins; ++i)
    EXPECT_EQ(2,l2.num_nnXsite(i));

  // Next-nearest neighbours
  EXPECT_EQ(36,l2.nnn.size());
  EXPECT_EQ(12,l2.nnnXsite.size());
  
  EXPECT_EQ(0,l2.nnn[0]);
  EXPECT_EQ(5,l2.nnn[1]);
  EXPECT_EQ(7,l2.nnn[2]);

  EXPECT_EQ(1,l2.nnn[3]);
  EXPECT_EQ(6,l2.nnn[4]);
  EXPECT_EQ(4,l2.nnn[5]);

  EXPECT_EQ(2,l2.nnn[6]);
  EXPECT_EQ(7,l2.nnn[7]);
  EXPECT_EQ(5,l2.nnn[8]);

  EXPECT_EQ(3,l2.nnn[9]);
  EXPECT_EQ(4,l2.nnn[10]);
  EXPECT_EQ(6,l2.nnn[11]);

  EXPECT_EQ(4,l2.nnn[12]);
  EXPECT_EQ(9,l2.nnn[13]);
  EXPECT_EQ(11,l2.nnn[14]);

  EXPECT_EQ(5,l2.nnn[15]);
  EXPECT_EQ(10,l2.nnn[16]);
  EXPECT_EQ(8,l2.nnn[17]);

  EXPECT_EQ(6,l2.nnn[18]);
  EXPECT_EQ(11,l2.nnn[19]);
  EXPECT_EQ(9,l2.nnn[20]);

  EXPECT_EQ(7,l2.nnn[21]);
  EXPECT_EQ(8,l2.nnn[22]);
  EXPECT_EQ(10,l2.nnn[23]);

  EXPECT_EQ(8,l2.nnn[24]);
  EXPECT_EQ(1,l2.nnn[25]);
  EXPECT_EQ(3,l2.nnn[26]);
  
  EXPECT_EQ(9,l2.nnn[27]);
  EXPECT_EQ(2,l2.nnn[28]);
  EXPECT_EQ(0,l2.nnn[29]);

  EXPECT_EQ(10,l2.nnn[30]);
  EXPECT_EQ(3,l2.nnn[31]);
  EXPECT_EQ(1,l2.nnn[32]);

  EXPECT_EQ(11,l2.nnn[33]);
  EXPECT_EQ(0,l2.nnn[34]);
  EXPECT_EQ(2,l2.nnn[35]);

  for (PetscInt i = 0; i < env.nspins; ++i)
    EXPECT_EQ(2,l2.num_nnnXsite(i));
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
