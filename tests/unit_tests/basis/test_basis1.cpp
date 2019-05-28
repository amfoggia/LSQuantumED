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

class BasisTest : public ::testing::Test {
protected:

  Basis A{env,0};
  
  std::vector<PetscInt> ref_int_basis1{7,11,13,14,19,
      21,22,25,26,28};
  std::vector<PetscInt> ref_int_basis2{35,37,38,41,42,
      44,49,50,52,56};
#ifdef DEVEL
  std::vector<boost::dynamic_bitset<>> ref_bit_basis1{boost::dynamic_bitset<>(6,7),
      boost::dynamic_bitset<>(6,11),
      boost::dynamic_bitset<>(6,13),
      boost::dynamic_bitset<>(6,14),
      boost::dynamic_bitset<>(6,19),
      boost::dynamic_bitset<>(6,21),
      boost::dynamic_bitset<>(6,22),
      boost::dynamic_bitset<>(6,25),
      boost::dynamic_bitset<>(6,26),
      boost::dynamic_bitset<>(6,28)};
  std::vector<boost::dynamic_bitset<>> ref_bit_basis2{boost::dynamic_bitset<>(6,35),
      boost::dynamic_bitset<>(6,37),
      boost::dynamic_bitset<>(6,38),
      boost::dynamic_bitset<>(6,41),
      boost::dynamic_bitset<>(6,42),
      boost::dynamic_bitset<>(6,44),
      boost::dynamic_bitset<>(6,49),
      boost::dynamic_bitset<>(6,50),
      boost::dynamic_bitset<>(6,52),
      boost::dynamic_bitset<>(6,56)};
#endif
};

TEST_F(BasisTest, defaultConstructor) {
  EXPECT_EQ(6, A.nspins);
  EXPECT_EQ(0, A.total_mag);
  EXPECT_EQ(3, A.nspins_up);
  EXPECT_EQ(20, A.size);  
  ASSERT_EQ(10, A.local_size);
  ASSERT_EQ(10, A.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
  ASSERT_EQ(10, A.bit_basis.size()); // meson.build is written to run with 2 processes
#endif
  
  if (A.mpi_rank == 0) {
    EXPECT_EQ(0, A.global_start_index);
    for (PetscInt i = 0; i < A.local_size; ++i) {
      EXPECT_EQ(ref_int_basis1.at(i), A.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis1.at(i), A.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
  else {
    EXPECT_EQ(10, A.global_start_index);
    for (PetscInt i = 0; i < A.local_size; ++i) {
      EXPECT_EQ(ref_int_basis2.at(i), A.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis2.at(i), A.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
}

TEST_F(BasisTest, copyConstructor) {
  Basis B{A};
  EXPECT_EQ(6, B.nspins);
  EXPECT_EQ(0, B.total_mag);
  EXPECT_EQ(3, B.nspins_up);
  EXPECT_EQ(20, B.size);  
  ASSERT_EQ(10, B.local_size);
  ASSERT_EQ(10, B.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
  ASSERT_EQ(10, B.bit_basis.size()); // meson.build is written to run with 2 processes
#endif

  if (B.mpi_rank == 0) {
    EXPECT_EQ(0, B.global_start_index);
      for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis1.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis1.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
      }
  }
  else {
    EXPECT_EQ(10, B.global_start_index);
    for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis2.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis2.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
}
  
TEST_F(BasisTest,moveConstructor) {
  Basis B{std::move(A)}; 
  EXPECT_EQ(6, B.nspins);
  EXPECT_EQ(0, B.total_mag);
  EXPECT_EQ(3, B.nspins_up);
  EXPECT_EQ(20, B.size);  
  ASSERT_EQ(10, B.local_size);
  ASSERT_EQ(10, B.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
  ASSERT_EQ(10, B.bit_basis.size()); // meson.build is written to run with 2 processes
#endif

  if (B.mpi_rank == 0) {
    EXPECT_EQ(0, B.global_start_index);
      for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis1.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis1.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
      }
  }
  else {
    EXPECT_EQ(10, B.global_start_index);
    for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis2.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis2.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
  
  EXPECT_EQ(0, A.int_basis.size());
  EXPECT_THROW(A.int_basis.at(1),std::out_of_range);
#ifdef DEVEL
  EXPECT_EQ(0, A.bit_basis.size());
  EXPECT_THROW(A.bit_basis.at(20),std::out_of_range);
#endif
}

TEST_F(BasisTest,copyAssignment) {
  Basis B{};
  B = A;
  EXPECT_EQ(6, B.nspins);
  EXPECT_EQ(0, B.total_mag);
  EXPECT_EQ(3, B.nspins_up);
  EXPECT_EQ(20, B.size);  
  ASSERT_EQ(10, B.local_size);
  ASSERT_EQ(10, B.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
  ASSERT_EQ(10, B.bit_basis.size()); // meson.build is written to run with 2 processes
#endif

  if (B.mpi_rank == 0) {
    EXPECT_EQ(0, B.global_start_index);
      for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis1.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis1.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
      }
  }
  else {
    EXPECT_EQ(10, B.global_start_index);
    for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis2.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis2.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }
}

TEST_F(BasisTest, moveAssignment) {
  Basis B{};
  B = std::move(A);
  EXPECT_EQ(6, B.nspins);
  EXPECT_EQ(0, B.total_mag);
  EXPECT_EQ(3, B.nspins_up);
  EXPECT_EQ(20, B.size);  
  ASSERT_EQ(10, B.local_size);
  ASSERT_EQ(10, B.int_basis.size()); // meson.build is written to run with 2 processes
#ifdef DEVEL
  ASSERT_EQ(10, B.bit_basis.size()); // meson.build is written to run with 2 processes
#endif

  if (B.mpi_rank == 0) {
    EXPECT_EQ(0, B.global_start_index);
      for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis1.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis1.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
      }
  }
  else {
    EXPECT_EQ(10, B.global_start_index);
    for (PetscInt i = 0; i < B.local_size; ++i) {
      EXPECT_EQ(ref_int_basis2.at(i), B.int_basis.at(i)) << "Vectors differ at index " << i;
#ifdef DEVEL
      EXPECT_EQ(ref_bit_basis2.at(i), B.bit_basis.at(i)) << "Vectors differ at index " << i;
#endif
    }
  }

  EXPECT_EQ(0, A.int_basis.size());
  EXPECT_THROW(A.int_basis.at(1),std::out_of_range);
#ifdef DEVEL
  EXPECT_EQ(0, A.bit_basis.size());
  EXPECT_THROW(A.bit_basis.at(20),std::out_of_range);
#endif
}

int main(int argc, char* argv[]) {
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
