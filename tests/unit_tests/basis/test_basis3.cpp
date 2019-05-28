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
  Basis AA{env,4};
  EXPECT_EQ(8, AA.nspins);
  EXPECT_EQ(4, AA.total_mag);
  EXPECT_EQ(8, AA.nspins_up);
  EXPECT_EQ(1, AA.size);
  EXPECT_EQ(1, AA.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, AA.mpi_rank);
  EXPECT_EQ(1, AA.local_size);
  EXPECT_EQ(0, AA.global_start_index);
  ASSERT_EQ(1, AA.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(1, AA.bit_basis.size());
#endif

  Basis BB{env,3};
  EXPECT_EQ(8, BB.nspins);
  EXPECT_EQ(3, BB.total_mag);
  EXPECT_EQ(7, BB.nspins_up);
  EXPECT_EQ(8, BB.size);
  EXPECT_EQ(1, BB.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, BB.mpi_rank);
  EXPECT_EQ(8, BB.local_size);
  EXPECT_EQ(0, BB.global_start_index); 
  ASSERT_EQ(8, BB.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(8, BB.bit_basis.size());
#endif

  Basis CC{env,2};
  EXPECT_EQ(8, CC.nspins);
  EXPECT_EQ(2, CC.total_mag);
  EXPECT_EQ(6, CC.nspins_up);
  EXPECT_EQ(28, CC.size);
  EXPECT_EQ(1, CC.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, CC.mpi_rank);
  EXPECT_EQ(28, CC.local_size);
  EXPECT_EQ(0, CC.global_start_index); 
  ASSERT_EQ(28, CC.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(28, CC.bit_basis.size());
#endif
  
  Basis DD{env,1};
  EXPECT_EQ(8, DD.nspins);
  EXPECT_EQ(1, DD.total_mag);
  EXPECT_EQ(5, DD.nspins_up);
  EXPECT_EQ(56, DD.size);
  EXPECT_EQ(1, DD.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, DD.mpi_rank);
  EXPECT_EQ(56, DD.local_size);
  EXPECT_EQ(0, DD.global_start_index); 
  ASSERT_EQ(56, DD.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(56, DD.bit_basis.size());
#endif
  
  Basis EE{env,-1};
  EXPECT_EQ(8, EE.nspins);
  EXPECT_EQ(-1, EE.total_mag);
  EXPECT_EQ(3, EE.nspins_up);
  EXPECT_EQ(56, EE.size);
  EXPECT_EQ(1, EE.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, EE.mpi_rank);
  EXPECT_EQ(56, EE.local_size);
  EXPECT_EQ(0, EE.global_start_index); 
  ASSERT_EQ(56, EE.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(56, EE.bit_basis.size());
#endif
  
  Basis FF{env,-2};
  EXPECT_EQ(8, FF.nspins);
  EXPECT_EQ(-2, FF.total_mag);
  EXPECT_EQ(2, FF.nspins_up);
  EXPECT_EQ(28, FF.size);
  EXPECT_EQ(1, FF.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, FF.mpi_rank);
  EXPECT_EQ(28, FF.local_size);
  EXPECT_EQ(0, FF.global_start_index); 
  ASSERT_EQ(28, FF.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(28, FF.bit_basis.size());
#endif
  
  Basis GG{env,-3};
  EXPECT_EQ(8, GG.nspins);
  EXPECT_EQ(-3, GG.total_mag);
  EXPECT_EQ(1, GG.nspins_up);
  EXPECT_EQ(8, GG.size);
  EXPECT_EQ(1, GG.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, GG.mpi_rank);
  EXPECT_EQ(8, GG.local_size);
  EXPECT_EQ(0, GG.global_start_index); 
  ASSERT_EQ(8, GG.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(8, GG.bit_basis.size());
#endif
  
  Basis HH{env,-4};
  EXPECT_EQ(8, HH.nspins);
  EXPECT_EQ(-4, HH.total_mag);
  EXPECT_EQ(0, HH.nspins_up);
  EXPECT_EQ(1, HH.size);
  EXPECT_EQ(1, HH.mpi_size); // meson.build with only one mpi process
  EXPECT_EQ(0, HH.mpi_rank);
  EXPECT_EQ(1, HH.local_size);
  EXPECT_EQ(0, HH.global_start_index); 
  ASSERT_EQ(1, HH.int_basis.size());
#ifdef DEVEL
  ASSERT_EQ(1, HH.bit_basis.size());
#endif
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
