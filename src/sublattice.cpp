#include "sublattice.hpp"

/* --------------------------------------------------------------------------- */

template<>
std::array<std::vector<PetscInt>,2> AF<2,chain1D>::construct_sublat() {

  std::array<std::vector<PetscInt>,2> sublat;
  
  for (PetscInt i = 0; i < nspins; ++i) 
    sublat[(i%2) ? 1 : 0].push_back(i);
  return sublat; 
}

/* --------------------------------------------------------------------------- */

template<>
std::array<std::vector<PetscInt>,2> AF<2,square2D>::construct_sublat() {

  std::array<std::vector<PetscInt>,2> sublat;
  std::vector<std::pair<PetscBool,short int>> tagged;
  short int current = 0, nn;
  PetscInt offset, spin;

  // Initialize tagged
  for (PetscInt i = 0; i < nspins; ++i)
    //    tagged.push_back(std::pair{PETSC_FALSE,0});
    tagged.push_back(std::pair<PetscBool,short int>{PETSC_FALSE,0});
  

  // Separate spins in sublattices
  for (PetscInt i = 0; i < nspins; ++i) {

    if (tagged[i].first == PETSC_FALSE) {
      tagged[i].first = PETSC_TRUE;
      tagged[i].second = current;
    }
    else
      current = tagged[i].second;

    nn = current ^ 1;

    offset = (this->lattice.num_neighbours() + 1) * i;
    for (PetscInt j = 0; j < this->lattice.num_neighbours(); ++j) {
      spin = this->lattice.nn[offset + j + 1];
      tagged[spin].first = PETSC_TRUE;
      tagged[spin].second = nn;
    }
  }

  // Fill in the sublattices
  for (PetscInt i = 0; i < nspins; ++ i) {
    if (tagged[i].second == 0)
      sublat[0].push_back(i);
    else
      sublat[1].push_back(i);
  } 
  return sublat;
}

/* --------------------------------------------------------------------------- */

template<>
std::array<std::vector<PetscInt>,4> Striped<4,square2D>::construct_sublat() {

  std::array<std::vector<PetscInt>,4> sublat;
  std::vector<std::pair<PetscBool,short int>> tagged;
  short int current;
  PetscInt neighA, neighB;
  PetscInt ix1,iy1;
  int current_spin;

  // First pair of sublattices ----------------------------------
  // Initialize tagged
  for (PetscInt i = 0; i < nspins; ++i)
    //    tagged.push_back(std::pair{PETSC_FALSE,0});
    tagged.push_back(std::pair<PetscBool,short int>{PETSC_FALSE,0});

  for (PetscInt ix = 0; ix < this->lattice.nspins_x; ++ix)
    for (PetscInt iy = 0; iy < this->lattice.nspins_y; ++iy) {
      current_spin = ix * this->lattice.nspins_y + iy;
      current = tagged[current_spin].second;

      iy1 = (iy+1)%this->lattice.nspins_y;
      neighA = ix * this->lattice.nspins_y + iy1;
      tagged[neighA].second = current ^ 1u;

      ix1 = (ix+1)%this->lattice.nspins_x;
      neighB = ix1 * this->lattice.nspins_y + iy;
      tagged[neighB].second = current;
    }

  // Fill in the sublattices
  for (PetscInt i = 0; i < nspins; ++ i) {
    if (tagged[i].second == 0)
      sublat[0].push_back(i);
    else
      sublat[1].push_back(i);
  }

  // Second pair of sublattices ----------------------------------
  // Initialize tagged
  for (PetscInt i = 0; i < nspins; ++i)
    //    tagged.push_back(std::pair{PETSC_FALSE,0});
    tagged.push_back(std::pair<PetscBool,short int>{PETSC_FALSE,0});

  for (PetscInt ix = 0; ix < this->lattice.nspins_x; ++ix)
    for (PetscInt iy = 0; iy < this->lattice.nspins_y; ++iy) {
      current_spin = ix * this->lattice.nspins_y + iy;
      current = tagged[current_spin].second;

      iy1 = (iy+1)%this->lattice.nspins_y;
      neighA = ix * this->lattice.nspins_y + iy1;
      tagged[neighA].second = current;

      ix1 = (ix+1)%this->lattice.nspins_x;
      neighB = ix1 * this->lattice.nspins_y + iy;
      tagged[neighB].second = current ^ 1u;
    }

  // Fill in the sublattices
  for (PetscInt i = 0; i < nspins; ++ i) {
    if (tagged[i].second == 0)
      sublat[2].push_back(i);
    else
      sublat[3].push_back(i);
  }
  return sublat;
}
