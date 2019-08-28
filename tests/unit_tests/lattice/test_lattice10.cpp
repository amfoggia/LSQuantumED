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

TEST(LatticeTest, honeycomb2D) {
  Environment env{_argc,_argv,16,help};
  
  honeycomb2D l1{env,4,4};

  EXPECT_EQ(lattice_type::honeycomb2D, l1.get_type());
  EXPECT_EQ(2,l1.num_nn());
  EXPECT_EQ(3,l1.num_nnn());

  // Nearest neighbours
  ASSERT_EQ(48,l1.nn.size());
  ASSERT_EQ(16,l1.nnXsite.size());

  EXPECT_EQ(0,l1.nn[0]);
  EXPECT_EQ(1,l1.nn[1]);
  EXPECT_EQ(4,l1.nn[2]);
  ASSERT_EQ(2,l1.num_nnXsite(0));

  EXPECT_EQ(1,l1.nn[3]);
  EXPECT_EQ(2,l1.nn[4]);
  EXPECT_EQ(-1,l1.nn[5]);
  ASSERT_EQ(1,l1.num_nnXsite(1));

  EXPECT_EQ(2,l1.nn[6]);
  EXPECT_EQ(3,l1.nn[7]);
  EXPECT_EQ(6,l1.nn[8]);
  ASSERT_EQ(2,l1.num_nnXsite(2));
  
  EXPECT_EQ(3,l1.nn[9]);
  EXPECT_EQ(0,l1.nn[10]);
  EXPECT_EQ(-1,l1.nn[11]);
  ASSERT_EQ(1,l1.num_nnXsite(3));

  EXPECT_EQ(4,l1.nn[12]);
  EXPECT_EQ(5,l1.nn[13]);
  EXPECT_EQ(-1,l1.nn[14]);
  ASSERT_EQ(1,l1.num_nnXsite(4));
  
  EXPECT_EQ(5,l1.nn[15]);
  EXPECT_EQ(6,l1.nn[16]);
  EXPECT_EQ(9,l1.nn[17]);
  ASSERT_EQ(2,l1.num_nnXsite(5));

  EXPECT_EQ(6,l1.nn[18]);
  EXPECT_EQ(7,l1.nn[19]);
  EXPECT_EQ(-1,l1.nn[20]);
  ASSERT_EQ(1,l1.num_nnXsite(6));

  EXPECT_EQ(7,l1.nn[21]);
  EXPECT_EQ(4,l1.nn[22]);
  EXPECT_EQ(11,l1.nn[23]);
  ASSERT_EQ(2,l1.num_nnXsite(7));

  EXPECT_EQ(8,l1.nn[24]);
  EXPECT_EQ(9,l1.nn[25]);
  EXPECT_EQ(12,l1.nn[26]);
  ASSERT_EQ(2,l1.num_nnXsite(8));
  
  EXPECT_EQ(9,l1.nn[27]);
  EXPECT_EQ(10,l1.nn[28]);
  EXPECT_EQ(-1,l1.nn[29]);
  ASSERT_EQ(1,l1.num_nnXsite(9));

  EXPECT_EQ(10,l1.nn[30]);
  EXPECT_EQ(11,l1.nn[31]);
  EXPECT_EQ(14,l1.nn[32]);
  ASSERT_EQ(2,l1.num_nnXsite(10));
  
  EXPECT_EQ(11,l1.nn[33]);
  EXPECT_EQ(8,l1.nn[34]);
  EXPECT_EQ(-1,l1.nn[35]);
  ASSERT_EQ(1,l1.num_nnXsite(11));

  EXPECT_EQ(12,l1.nn[36]);
  EXPECT_EQ(13,l1.nn[37]);
  EXPECT_EQ(-1,l1.nn[38]);
  ASSERT_EQ(1,l1.num_nnXsite(12));

  EXPECT_EQ(13,l1.nn[39]);
  EXPECT_EQ(14,l1.nn[40]);
  EXPECT_EQ(1,l1.nn[41]);
  ASSERT_EQ(2,l1.num_nnXsite(13));

  EXPECT_EQ(14,l1.nn[42]);
  EXPECT_EQ(15,l1.nn[43]);
  EXPECT_EQ(-1,l1.nn[44]);
  ASSERT_EQ(1,l1.num_nnXsite(14));
  
  EXPECT_EQ(15,l1.nn[45]);
  EXPECT_EQ(12,l1.nn[46]);
  EXPECT_EQ(3,l1.nn[47]);
  ASSERT_EQ(2,l1.num_nnXsite(15));

  // Next-nearest neighbours
  ASSERT_EQ(64,l1.nnn.size());
  ASSERT_EQ(16,l1.nnnXsite.size());

  EXPECT_EQ(0,l1.nnn[0]);
  EXPECT_EQ(2,l1.nnn[1]);
  EXPECT_EQ(5,l1.nnn[2]);
  EXPECT_EQ(7,l1.nnn[3]);
  
  EXPECT_EQ(1,l1.nnn[4]);
  EXPECT_EQ(3,l1.nnn[5]);
  EXPECT_EQ(6,l1.nnn[6]);
  EXPECT_EQ(4,l1.nnn[7]);
  
  EXPECT_EQ(2,l1.nnn[8]);
  EXPECT_EQ(0,l1.nnn[9]);
  EXPECT_EQ(7,l1.nnn[10]);
  EXPECT_EQ(5,l1.nnn[11]);

  EXPECT_EQ(3,l1.nnn[12]);
  EXPECT_EQ(1,l1.nnn[13]);
  EXPECT_EQ(4,l1.nnn[14]);
  EXPECT_EQ(6,l1.nnn[15]);
  
  EXPECT_EQ(4,l1.nnn[16]);
  EXPECT_EQ(6,l1.nnn[17]);
  EXPECT_EQ(9,l1.nnn[18]);
  EXPECT_EQ(11,l1.nnn[19]);
  
  EXPECT_EQ(5,l1.nnn[20]);
  EXPECT_EQ(7,l1.nnn[21]);
  EXPECT_EQ(10,l1.nnn[22]);
  EXPECT_EQ(8,l1.nnn[23]);

  EXPECT_EQ(6,l1.nnn[24]);
  EXPECT_EQ(4,l1.nnn[25]);
  EXPECT_EQ(11,l1.nnn[26]);
  EXPECT_EQ(9,l1.nnn[27]);
  
  EXPECT_EQ(7,l1.nnn[28]);
  EXPECT_EQ(5,l1.nnn[29]);
  EXPECT_EQ(8,l1.nnn[30]);
  EXPECT_EQ(10,l1.nnn[31]);
  
  EXPECT_EQ(8,l1.nnn[32]);
  EXPECT_EQ(10,l1.nnn[33]);
  EXPECT_EQ(13,l1.nnn[34]);
  EXPECT_EQ(15,l1.nnn[35]);

  EXPECT_EQ(9,l1.nnn[36]);
  EXPECT_EQ(11,l1.nnn[37]);
  EXPECT_EQ(14,l1.nnn[38]);
  EXPECT_EQ(12,l1.nnn[39]);
  
  EXPECT_EQ(10,l1.nnn[40]);
  EXPECT_EQ(8,l1.nnn[41]);
  EXPECT_EQ(15,l1.nnn[42]);
  EXPECT_EQ(13,l1.nnn[43]);
  
  EXPECT_EQ(11,l1.nnn[44]);
  EXPECT_EQ(9,l1.nnn[45]);
  EXPECT_EQ(12,l1.nnn[46]);
  EXPECT_EQ(14,l1.nnn[47]);

  EXPECT_EQ(12,l1.nnn[48]);
  EXPECT_EQ(14,l1.nnn[49]);
  EXPECT_EQ(1,l1.nnn[50]);
  EXPECT_EQ(3,l1.nnn[51]);
  
  EXPECT_EQ(13,l1.nnn[52]);
  EXPECT_EQ(15,l1.nnn[53]);
  EXPECT_EQ(2,l1.nnn[54]);
  EXPECT_EQ(0,l1.nnn[55]);
  
  EXPECT_EQ(14,l1.nnn[56]);
  EXPECT_EQ(12,l1.nnn[57]);
  EXPECT_EQ(3,l1.nnn[58]);
  EXPECT_EQ(1,l1.nnn[59]);

  EXPECT_EQ(15,l1.nnn[60]);
  EXPECT_EQ(13,l1.nnn[61]);
  EXPECT_EQ(0,l1.nnn[62]);  
  EXPECT_EQ(2,l1.nnn[63]);

  for (PetscInt i = 0; i < 16; ++i)
    ASSERT_EQ(3,l1.num_nnnXsite(i));
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
