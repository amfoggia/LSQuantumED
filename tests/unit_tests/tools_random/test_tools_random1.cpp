#include "random_disorder.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

int _argc;
char ** _argv;

class HamiltonianTestEnv : public ::testing::Environment {
protected:
  virtual void TearDown() {
      ::testing::internal::CaptureStdout();
      EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};

TEST(ToolsRandom, Constructor) {
  {
    Tools::RandomDisorder rd{6, PETSC_TRUE, -0.005, 0.005, 10};
    EXPECT_EQ(6, rd.nspins);
    EXPECT_NEAR(-0.005, rd.min, 1e-13);
    EXPECT_NEAR(0.005, rd.max, 1e-13);
    EXPECT_EQ(10, rd.rep);

    EXPECT_EQ(134, rd.srand_seq[0]);
    EXPECT_EQ(136, rd.srand_seq[1]);
    EXPECT_EQ(451, rd.srand_seq[2]);
    EXPECT_EQ(21, rd.srand_seq[3]);
    EXPECT_EQ(351, rd.srand_seq[4]);
    EXPECT_EQ(912, rd.srand_seq[5]);
    EXPECT_EQ(471, rd.srand_seq[6]);
    EXPECT_EQ(74, rd.srand_seq[7]);
    EXPECT_EQ(570, rd.srand_seq[8]);
    EXPECT_EQ(635, rd.srand_seq[9]);

    rd.hi_init(2);
    EXPECT_NEAR(0.00033296339931027097, rd.hi[0] , 1e-13);
    EXPECT_NEAR(-0.0024819526282803856, rd.hi[1] , 1e-13);
    EXPECT_NEAR(-0.0043726977344249567, rd.hi[2] , 1e-13);
    EXPECT_NEAR(-0.0043282870609637177, rd.hi[3] , 1e-13);
    EXPECT_NEAR(-0.0018033511217890774, rd.hi[4] , 1e-13);
    EXPECT_NEAR(-0.0002136219343576908, rd.hi[5] , 1e-13);
  }
}

int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}
