#include "basis.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Basis class\n\n";
int _argc;
char ** _argv;

Environment env{_argc,_argv,8,help};

TEST(MembersTest, nspins) {
  ASSERT_THROW(Basis(env,4), std::invalid_argument);
  
  
  Basis BB{env,3};
  EXPECT_EQ(8, BB.nspins);
  EXPECT_EQ(3, BB.total_mag);
  EXPECT_EQ(7, BB.nspins_up);
  EXPECT_EQ(8, BB.size);
  EXPECT_EQ(5, BB.mpi_size); // meson.build with only one mpi process

  if (BB.mpi_rank == 0) {
    EXPECT_EQ(2, BB.local_size);
    EXPECT_EQ(0, BB.global_start_index);
    ASSERT_EQ(2, BB.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(2, BB.bit_basis.size());
#endif
  }
  else if (BB.mpi_rank == 1) {
    EXPECT_EQ(2, BB.local_size);
    EXPECT_EQ(2, BB.global_start_index);
    ASSERT_EQ(2, BB.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(2, BB.bit_basis.size());
#endif
  }
  else if (BB.mpi_rank == 2) {
    EXPECT_EQ(2, BB.local_size);
    EXPECT_EQ(4, BB.global_start_index);
    ASSERT_EQ(2, BB.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(2, BB.bit_basis.size());
#endif
  }
  else if (BB.mpi_rank == 3) {
    EXPECT_EQ(1, BB.local_size);
    EXPECT_EQ(6, BB.global_start_index);
    ASSERT_EQ(1, BB.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(1, BB.bit_basis.size());
#endif
  }
  else {
    EXPECT_EQ(1, BB.local_size);
    EXPECT_EQ(7, BB.global_start_index);
    ASSERT_EQ(1, BB.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(1, BB.bit_basis.size());
#endif
  }
  
  Basis CC{env,2};
  EXPECT_EQ(8, CC.nspins);
  EXPECT_EQ(2, CC.total_mag);
  EXPECT_EQ(6, CC.nspins_up);
  EXPECT_EQ(28, CC.size);
  EXPECT_EQ(5, CC.mpi_size); // meson.build with only one mpi process
  
  if (CC.mpi_rank == 0) {
    EXPECT_EQ(6, CC.local_size);
    EXPECT_EQ(0, CC.global_start_index);
    ASSERT_EQ(6, CC.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(6, CC.bit_basis.size());
#endif
  }
  else if (CC.mpi_rank == 1) {
    EXPECT_EQ(6, CC.local_size);
    EXPECT_EQ(6, CC.global_start_index);
    ASSERT_EQ(6, CC.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(6, CC.bit_basis.size());
#endif
  }
  else if (CC.mpi_rank == 2) {
    EXPECT_EQ(6, CC.local_size);
    EXPECT_EQ(12, CC.global_start_index);
    ASSERT_EQ(6, CC.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(6, CC.bit_basis.size());
#endif
  }
  else if (CC.mpi_rank == 3) {
    EXPECT_EQ(5, CC.local_size);
    EXPECT_EQ(18, CC.global_start_index);
    ASSERT_EQ(5, CC.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(5, CC.bit_basis.size());
#endif
  }
  else {
    EXPECT_EQ(5, CC.local_size);
    EXPECT_EQ(23, CC.global_start_index);
    ASSERT_EQ(5, CC.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(5, CC.bit_basis.size());
#endif
  }

  Basis DD{env,1};
  EXPECT_EQ(8, DD.nspins);
  EXPECT_EQ(1, DD.total_mag);
  EXPECT_EQ(5, DD.nspins_up);
  EXPECT_EQ(56, DD.size);
  EXPECT_EQ(5, DD.mpi_size); // meson.build with only one mpi process
  
  if (DD.mpi_rank == 0) {
    EXPECT_EQ(12, DD.local_size);
    EXPECT_EQ(0, DD.global_start_index);
    ASSERT_EQ(12, DD.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(12, DD.bit_basis.size());
#endif
  }
  else if (DD.mpi_rank == 1) {
    EXPECT_EQ(11, DD.local_size);
    EXPECT_EQ(12, DD.global_start_index);
    ASSERT_EQ(11, DD.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(11, DD.bit_basis.size());
#endif
  }
  else if (DD.mpi_rank == 2) {
    EXPECT_EQ(11, DD.local_size);
    EXPECT_EQ(23, DD.global_start_index);
    ASSERT_EQ(11, DD.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(11, DD.bit_basis.size());
#endif
  }
  else if (DD.mpi_rank == 3) {
    EXPECT_EQ(11, DD.local_size);
    EXPECT_EQ(34, DD.global_start_index);
    ASSERT_EQ(11, DD.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(11, DD.bit_basis.size());
#endif
  }
  else {
    EXPECT_EQ(11, DD.local_size);
    EXPECT_EQ(45, DD.global_start_index);
    ASSERT_EQ(11, DD.int_basis.size());
#ifdef DEVEL
    ASSERT_EQ(11, DD.bit_basis.size());
#endif
  }
}

TEST(MembersTest, totalMag) {
  EXPECT_THROW(Basis a(env,5), std::invalid_argument);
  EXPECT_THROW(Basis a(env,-6), std::invalid_argument);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
