#include "environment.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Environment class\n\n";
int _argc;
char ** _argv;

TEST(EnvironmentTest, ConstrDestrManySpins) {
  ASSERT_THROW(Environment testEnv(_argc, _argv, 37, help), std::invalid_argument);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
