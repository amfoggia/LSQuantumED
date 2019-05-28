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

Environment env{_argc,_argv,6,help};

TEST(BasisTest, BasisConstructionAllCases61) {
  Basis K_61{env,1};
  std::vector<PetscInt> ref_int_basis_61_1{15,23,27,29,30,
      39,43,45};
  std::vector<PetscInt> ref_int_basis_61_2{46,51,53,54,57,58,60};
#ifdef DEVEL
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_61_1{boost::dynamic_bitset<>(6,15),
      boost::dynamic_bitset<>(6,23),
      boost::dynamic_bitset<>(6,27),
      boost::dynamic_bitset<>(6,29),
      boost::dynamic_bitset<>(6,30),
      boost::dynamic_bitset<>(6,39),
      boost::dynamic_bitset<>(6,43),
      boost::dynamic_bitset<>(6,45)};
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_61_2{boost::dynamic_bitset<>(6,46),
      boost::dynamic_bitset<>(6,51),
      boost::dynamic_bitset<>(6,53),
      boost::dynamic_bitset<>(6,54),
      boost::dynamic_bitset<>(6,57),
      boost::dynamic_bitset<>(6,58),
      boost::dynamic_bitset<>(6,60)};
#endif
  EXPECT_EQ(15, K_61.size);
  
  if (K_61.mpi_rank == 0) {
    ASSERT_EQ(8, K_61.local_size);
    EXPECT_EQ(0, K_61.global_start_index);
    ASSERT_EQ(8, K_61.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(8, K_61.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_61.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_61_1.at(i), K_61.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_61_1.at(i), K_61.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
  else {
    ASSERT_EQ(7, K_61.local_size);
    EXPECT_EQ(8, K_61.global_start_index);
    ASSERT_EQ(7, K_61.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(7, K_61.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_61.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_61_2.at(i), K_61.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_61_2.at(i), K_61.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
}

TEST(BasisTest, BasisConstructionAllCases6m1) {
  Basis K_6m1{env,-1};
  std::vector<PetscInt> ref_int_basis_6m1_1{3,5,6,9,10,
      12,17,18};
  std::vector<PetscInt> ref_int_basis_6m1_2{20,24,
      33,34,36,40,48};
#ifdef DEVEL
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_6m1_1{boost::dynamic_bitset<>(6,3),
      boost::dynamic_bitset<>(6,5),
      boost::dynamic_bitset<>(6,6),
      boost::dynamic_bitset<>(6,9),
      boost::dynamic_bitset<>(6,10),
      boost::dynamic_bitset<>(6,12),
      boost::dynamic_bitset<>(6,17),
      boost::dynamic_bitset<>(6,18)};
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_6m1_2{boost::dynamic_bitset<>(6,20),
      boost::dynamic_bitset<>(6,24),
      boost::dynamic_bitset<>(6,33),
      boost::dynamic_bitset<>(6,34),
      boost::dynamic_bitset<>(6,36),
      boost::dynamic_bitset<>(6,40),
      boost::dynamic_bitset<>(6,48)};
#endif
  
  EXPECT_EQ(15, K_6m1.size);  

  if (K_6m1.mpi_rank == 0) {
    ASSERT_EQ(8, K_6m1.local_size);
    EXPECT_EQ(0, K_6m1.global_start_index);
    ASSERT_EQ(8, K_6m1.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(8, K_6m1.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_6m1.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_6m1_1.at(i), K_6m1.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_6m1_1.at(i), K_6m1.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
  else {
    ASSERT_EQ(7, K_6m1.local_size);
    EXPECT_EQ(8, K_6m1.global_start_index);
    ASSERT_EQ(7, K_6m1.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(7, K_6m1.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_6m1.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_6m1_2.at(i), K_6m1.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_6m1_2.at(i), K_6m1.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
}

TEST(BasisTest, BasisConstructionAllCases62) {
  Basis K_62{env,2};
  std::vector<PetscInt> ref_int_basis_62_1{31,47,55};
  std::vector<PetscInt> ref_int_basis_62_2{59,61,62};
#ifdef DEVEL
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_62_1{boost::dynamic_bitset<>(6,31),
      boost::dynamic_bitset<>(6,47),
      boost::dynamic_bitset<>(6,55)};
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_62_2{boost::dynamic_bitset<>(6,59),
      boost::dynamic_bitset<>(6,61),
      boost::dynamic_bitset<>(6,62)};
#endif
  
  EXPECT_EQ(6, K_62.size);  

  if (K_62.mpi_rank == 0) {
    ASSERT_EQ(3, K_62.local_size);
    EXPECT_EQ(0, K_62.global_start_index);
    ASSERT_EQ(3, K_62.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(3, K_62.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_62.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_62_1.at(i), K_62.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_62_1.at(i), K_62.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
  else {
    ASSERT_EQ(3, K_62.local_size);
    EXPECT_EQ(3, K_62.global_start_index);
    ASSERT_EQ(3, K_62.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(3, K_62.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_62.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_62_2.at(i), K_62.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_62_2.at(i), K_62.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
}

TEST(BasisTest, BasisConstructionAllCases6m2) {
  Basis K_6m2{env,-2};
  std::vector<PetscInt> ref_int_basis_6m2_1{1,2,4};
  std::vector<PetscInt> ref_int_basis_6m2_2{8,16,32};
#ifdef DEVEL
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_6m2_1{boost::dynamic_bitset<>(6,1),
      boost::dynamic_bitset<>(6,2),
      boost::dynamic_bitset<>(6,4)};
  std::vector<boost::dynamic_bitset<>> ref_bit_basis_6m2_2{boost::dynamic_bitset<>(6,8),
      boost::dynamic_bitset<>(6,16),
      boost::dynamic_bitset<>(6,32)};
#endif
  
  EXPECT_EQ(6, K_6m2.size);  

  if (K_6m2.mpi_rank == 0) {
    ASSERT_EQ(3, K_6m2.local_size);
    EXPECT_EQ(0, K_6m2.global_start_index);
    ASSERT_EQ(3, K_6m2.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(3, K_6m2.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_6m2.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_6m2_1.at(i), K_6m2.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_6m2_1.at(i), K_6m2.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
  else {
    ASSERT_EQ(3, K_6m2.local_size);
    EXPECT_EQ(3, K_6m2.global_start_index);
    ASSERT_EQ(3, K_6m2.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
    ASSERT_EQ(3, K_6m2.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
    for (PetscInt i = 0; i < K_6m2.local_size; ++i) {
      EXPECT_EQ(ref_int_basis_6m2_2.at(i), K_6m2.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis_6m2_2.at(i), K_6m2.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
}

TEST(BasisTest, BasisConstructionAllCases63) {
  ASSERT_THROW(Basis(env,3),std::invalid_argument);
}

TEST(BasisTest, BasisConstructionAllCases6m3) {
  ASSERT_THROW(Basis(env,-3),std::invalid_argument);
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
