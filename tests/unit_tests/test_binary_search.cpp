#include <gtest/gtest.h>
#include "hamiltonian.hpp"

TEST(HamiltonianTest, BinarySearch) {
  std::vector<PetscInt> array;
  int x;

  std::srand(7);

  for (int i = 0; i < 30; ++i) {
    x = (std::rand() % static_cast<int>(30 + 1));
    array.push_back(x);
  }

  std::sort(array.begin(), array.end());

  EXPECT_EQ(12, HamiltHelper::binary_search(array, 0, 29, 15));
  EXPECT_EQ(1, HamiltHelper::binary_search(array, 0, 29, 4));
  EXPECT_EQ(24, HamiltHelper::binary_search(array, 0, 29, 26));
  
}


