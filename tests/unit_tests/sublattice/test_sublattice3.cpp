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

TEST(SublatticeTest, GetSublatsHoneycomb2D) {
  Environment env{_argc,_argv,16,help};

/* --------------------------------------------------------------------------- */
  
  honeycomb2D l1{env,4,4};
  AF<honeycomb2D> sublat1{l1,2};
  EXPECT_EQ(2,sublat1.get_size());
  
  EXPECT_EQ(0,sublat1.get_sl(0)[0]);
  EXPECT_EQ(2,sublat1.get_sl(0)[1]);
  EXPECT_EQ(5,sublat1.get_sl(0)[2]);
  EXPECT_EQ(7,sublat1.get_sl(0)[3]);
  EXPECT_EQ(8,sublat1.get_sl(0)[4]);
  EXPECT_EQ(10,sublat1.get_sl(0)[5]);
  EXPECT_EQ(13,sublat1.get_sl(0)[6]);
  EXPECT_EQ(15,sublat1.get_sl(0)[7]);
  
  EXPECT_EQ(1,sublat1.get_sl(1)[0]);
  EXPECT_EQ(3,sublat1.get_sl(1)[1]);
  EXPECT_EQ(4,sublat1.get_sl(1)[2]);
  EXPECT_EQ(6,sublat1.get_sl(1)[3]);
  EXPECT_EQ(9,sublat1.get_sl(1)[4]);
  EXPECT_EQ(11,sublat1.get_sl(1)[5]);
  EXPECT_EQ(12,sublat1.get_sl(1)[6]);
  EXPECT_EQ(14,sublat1.get_sl(1)[7]);

  Striped<honeycomb2D> sublat11{l1,6};
  EXPECT_EQ(6,sublat11.get_size());
  
  EXPECT_EQ(0,sublat11.get_sl(0)[0]);
  EXPECT_EQ(2,sublat11.get_sl(0)[1]);
  EXPECT_EQ(4,sublat11.get_sl(0)[2]);
  EXPECT_EQ(6,sublat11.get_sl(0)[3]);
  EXPECT_EQ(8,sublat11.get_sl(0)[4]);
  EXPECT_EQ(10,sublat11.get_sl(0)[5]);
  EXPECT_EQ(12,sublat11.get_sl(0)[6]);
  EXPECT_EQ(14,sublat11.get_sl(0)[7]);
  
  EXPECT_EQ(1,sublat11.get_sl(1)[0]);
  EXPECT_EQ(3,sublat11.get_sl(1)[1]);
  EXPECT_EQ(5,sublat11.get_sl(1)[2]);
  EXPECT_EQ(7,sublat11.get_sl(1)[3]);
  EXPECT_EQ(9,sublat11.get_sl(1)[4]);
  EXPECT_EQ(11,sublat11.get_sl(1)[5]);
  EXPECT_EQ(13,sublat11.get_sl(1)[6]);
  EXPECT_EQ(15,sublat11.get_sl(1)[7]);

  // EXPECT_EQ(0,sublat11.get_sl(2)[0]);
  // EXPECT_EQ(1,sublat11.get_sl(2)[1]);
  // EXPECT_EQ(2,sublat11.get_sl(2)[2]);
  // EXPECT_EQ(3,sublat11.get_sl(2)[3]);
  // EXPECT_EQ(8,sublat11.get_sl(2)[4]);
  // EXPECT_EQ(9,sublat11.get_sl(2)[5]);
  // EXPECT_EQ(10,sublat11.get_sl(2)[6]);
  // EXPECT_EQ(11,sublat11.get_sl(2)[7]);

  // EXPECT_EQ(4,sublat11.get_sl(3)[0]);
  // EXPECT_EQ(5,sublat11.get_sl(3)[1]);
  // EXPECT_EQ(6,sublat11.get_sl(3)[2]);
  // EXPECT_EQ(7,sublat11.get_sl(3)[3]);
  // EXPECT_EQ(12,sublat11.get_sl(3)[4]);
  // EXPECT_EQ(13,sublat11.get_sl(3)[5]);
  // EXPECT_EQ(14,sublat11.get_sl(3)[6]);
  // EXPECT_EQ(15,sublat11.get_sl(3)[7]);

/* --------------------------------------------------------------------------- */
  
  honeycomb2D l2{env,2,8};
  AF<honeycomb2D> sublat2{l2,4};
  EXPECT_EQ(2,sublat2.get_size());
  
  EXPECT_EQ(0,sublat2.get_sl(0)[0]);
  EXPECT_EQ(2,sublat2.get_sl(0)[1]);
  EXPECT_EQ(4,sublat2.get_sl(0)[2]);
  EXPECT_EQ(6,sublat2.get_sl(0)[3]);
  EXPECT_EQ(9,sublat2.get_sl(0)[4]);
  EXPECT_EQ(11,sublat2.get_sl(0)[5]);
  EXPECT_EQ(13,sublat2.get_sl(0)[6]);
  EXPECT_EQ(15,sublat2.get_sl(0)[7]);
 
  EXPECT_EQ(1,sublat2.get_sl(1)[0]);
  EXPECT_EQ(3,sublat2.get_sl(1)[1]);
  EXPECT_EQ(5,sublat2.get_sl(1)[2]);
  EXPECT_EQ(7,sublat2.get_sl(1)[3]);
  EXPECT_EQ(8,sublat2.get_sl(1)[4]);
  EXPECT_EQ(10,sublat2.get_sl(1)[5]);
  EXPECT_EQ(12,sublat2.get_sl(1)[6]);
  EXPECT_EQ(14,sublat2.get_sl(1)[7]);

  // Striped<honeycomb2D> sublat22{l2};
  // EXPECT_EQ(4,sublat22.get_size());
  
  // EXPECT_EQ(0,sublat22.get_sl(0)[0]);
  // EXPECT_EQ(2,sublat22.get_sl(0)[1]);
  // EXPECT_EQ(4,sublat22.get_sl(0)[2]);
  // EXPECT_EQ(6,sublat22.get_sl(0)[3]);
  // EXPECT_EQ(8,sublat22.get_sl(0)[4]);
  // EXPECT_EQ(10,sublat22.get_sl(0)[5]);
  // EXPECT_EQ(12,sublat22.get_sl(0)[6]);
  // EXPECT_EQ(14,sublat22.get_sl(0)[7]);

  // EXPECT_EQ(1,sublat22.get_sl(1)[0]);
  // EXPECT_EQ(3,sublat22.get_sl(1)[1]);
  // EXPECT_EQ(5,sublat22.get_sl(1)[2]);
  // EXPECT_EQ(7,sublat22.get_sl(1)[3]);
  // EXPECT_EQ(9,sublat22.get_sl(1)[4]);
  // EXPECT_EQ(11,sublat22.get_sl(1)[5]);
  // EXPECT_EQ(13,sublat22.get_sl(1)[6]);
  // EXPECT_EQ(15,sublat22.get_sl(1)[7]);

  // EXPECT_EQ(0,sublat22.get_sl(2)[0]);
  // EXPECT_EQ(1,sublat22.get_sl(2)[1]);
  // EXPECT_EQ(2,sublat22.get_sl(2)[2]);
  // EXPECT_EQ(3,sublat22.get_sl(2)[3]);
  // EXPECT_EQ(4,sublat22.get_sl(2)[4]);
  // EXPECT_EQ(5,sublat22.get_sl(2)[5]);
  // EXPECT_EQ(6,sublat22.get_sl(2)[6]);
  // EXPECT_EQ(7,sublat22.get_sl(2)[7]);

  // EXPECT_EQ(8,sublat22.get_sl(3)[0]);
  // EXPECT_EQ(9,sublat22.get_sl(3)[1]);
  // EXPECT_EQ(10,sublat22.get_sl(3)[2]);
  // EXPECT_EQ(11,sublat22.get_sl(3)[3]);
  // EXPECT_EQ(12,sublat22.get_sl(3)[4]);
  // EXPECT_EQ(13,sublat22.get_sl(3)[5]);
  // EXPECT_EQ(14,sublat22.get_sl(3)[6]);
  // EXPECT_EQ(15,sublat22.get_sl(3)[7]);

/* --------------------------------------------------------------------------- */
  
  honeycomb2D l3{env,8,2};
  AF<honeycomb2D> sublat3{l3,4};
  EXPECT_EQ(2,sublat3.get_size());
  
  EXPECT_EQ(0,sublat3.get_sl(0)[0]);
  EXPECT_EQ(3,sublat3.get_sl(0)[1]);
  EXPECT_EQ(4,sublat3.get_sl(0)[2]);
  EXPECT_EQ(7,sublat3.get_sl(0)[3]);
  EXPECT_EQ(8,sublat3.get_sl(0)[4]);
  EXPECT_EQ(11,sublat3.get_sl(0)[5]);
  EXPECT_EQ(12,sublat3.get_sl(0)[6]);
  EXPECT_EQ(15,sublat3.get_sl(0)[7]);

  EXPECT_EQ(1,sublat3.get_sl(1)[0]);
  EXPECT_EQ(2,sublat3.get_sl(1)[1]);
  EXPECT_EQ(5,sublat3.get_sl(1)[2]);
  EXPECT_EQ(6,sublat3.get_sl(1)[3]);
  EXPECT_EQ(9,sublat3.get_sl(1)[4]);
  EXPECT_EQ(10,sublat3.get_sl(1)[5]);
  EXPECT_EQ(13,sublat3.get_sl(1)[6]);
  EXPECT_EQ(14,sublat3.get_sl(1)[7]);

  // Striped<honeycomb2D> sublat33{l3};
  // EXPECT_EQ(4,sublat33.get_size());
  
  // EXPECT_EQ(0,sublat33.get_sl(0)[0]);
  // EXPECT_EQ(2,sublat33.get_sl(0)[1]);
  // EXPECT_EQ(4,sublat33.get_sl(0)[2]);
  // EXPECT_EQ(6,sublat33.get_sl(0)[3]);
  // EXPECT_EQ(8,sublat33.get_sl(0)[4]);
  // EXPECT_EQ(10,sublat33.get_sl(0)[5]);
  // EXPECT_EQ(12,sublat33.get_sl(0)[6]);
  // EXPECT_EQ(14,sublat33.get_sl(0)[7]);

  // EXPECT_EQ(1,sublat33.get_sl(1)[0]);
  // EXPECT_EQ(3,sublat33.get_sl(1)[1]);
  // EXPECT_EQ(5,sublat33.get_sl(1)[2]);
  // EXPECT_EQ(7,sublat33.get_sl(1)[3]);
  // EXPECT_EQ(9,sublat33.get_sl(1)[4]);
  // EXPECT_EQ(11,sublat33.get_sl(1)[5]);
  // EXPECT_EQ(13,sublat33.get_sl(1)[6]);
  // EXPECT_EQ(15,sublat33.get_sl(1)[7]);

  // EXPECT_EQ(0,sublat33.get_sl(2)[0]);
  // EXPECT_EQ(1,sublat33.get_sl(2)[1]);
  // EXPECT_EQ(4,sublat33.get_sl(2)[2]);
  // EXPECT_EQ(5,sublat33.get_sl(2)[3]);
  // EXPECT_EQ(8,sublat33.get_sl(2)[4]);
  // EXPECT_EQ(9,sublat33.get_sl(2)[5]);
  // EXPECT_EQ(12,sublat33.get_sl(2)[6]);
  // EXPECT_EQ(13,sublat33.get_sl(2)[7]);

  // EXPECT_EQ(2,sublat33.get_sl(3)[0]);
  // EXPECT_EQ(3,sublat33.get_sl(3)[1]);
  // EXPECT_EQ(6,sublat33.get_sl(3)[2]);
  // EXPECT_EQ(7,sublat33.get_sl(3)[3]);
  // EXPECT_EQ(10,sublat33.get_sl(3)[4]);
  // EXPECT_EQ(11,sublat33.get_sl(3)[5]);
  // EXPECT_EQ(14,sublat33.get_sl(3)[6]);
  // EXPECT_EQ(15,sublat33.get_sl(3)[7]);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
