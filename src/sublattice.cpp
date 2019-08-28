#include "sublattice.hpp"

/* --------------------------------------------------------------------------- */

template<>
std::vector<std::vector<PetscInt>> AF<chain1D>::construct_sublat() {

  std::vector<std::vector<PetscInt>> sublat;

  for (PetscInt s = 0; s < size; ++s)
    sublat.push_back(std::vector<PetscInt>{});
  
  for (PetscInt i = 0; i < nspins; ++i) 
    sublat[(i%2) ? 1 : 0].push_back(i);
  return sublat; 
}

template<>
AF<chain1D>::AF(chain1D& m_lattice,
		PetscInt m_size)
  : Sublattice<chain1D>{m_lattice},
  nspins{m_lattice.num_spins()},
  size{2}
  {sublat = construct_sublat();}

/* --------------------------------------------------------------------------- */

template<>
std::vector<std::vector<PetscInt>> AF<square2D>::construct_sublat() {

  std::vector<std::vector<PetscInt>> sublat;
  short int tagx, tagy;
  PetscInt current;
  
  for (PetscInt s = 0; s < size; ++s)
    sublat.push_back(std::vector<PetscInt>{});

  tagx = 0;
  for (PetscInt ix = 0; ix < this->lattice.nspins_x; ++ix) {
    tagy = tagx;
    for (PetscInt iy = 0; iy < this->lattice.nspins_y; ++iy) {
      current = ix * this->lattice.nspins_y + iy;
      sublat[(tagy%2) ? 1 : 0].push_back(current);
      tagy ^= 1;
    }
    tagx ^= 1;
  }
  return sublat;
}

template<>
AF<square2D>::AF(square2D& m_lattice,
		 PetscInt m_size)
  : Sublattice<square2D>{m_lattice},
  nspins{m_lattice.num_spins()},
  size{2}
  {sublat = construct_sublat();}

/* --------------------------------------------------------------------------- */

template<>
std::vector<std::vector<PetscInt>> Striped<square2D>::construct_sublat() {

  std::vector<std::vector<PetscInt>> sublat;
  PetscInt current;
  short int tagx, tagy;
  
  for (PetscInt s = 0; s < size; ++s)
    sublat.push_back(std::vector<PetscInt>{});

  tagx = 0;
  for (PetscInt ix = 0; ix < this->lattice.nspins_x; ++ix) {
    tagy = 0;
    for (PetscInt iy = 0; iy < this->lattice.nspins_y; ++iy) {
      current = ix * this->lattice.nspins_y + iy;
      sublat[(tagy%2) ? 1 : 0].push_back(current);
      sublat[(tagx%2) ? 3 : 2].push_back(current);
      tagy ^= 1;
    }
    tagx ^= 1;
  }
  return sublat;
}

template<>
Striped<square2D>::Striped(square2D& m_lattice,
			   PetscInt m_size)
  : Sublattice<square2D>{m_lattice},
  nspins{m_lattice.nspins},
  size{4}
{sublat = construct_sublat();}

/* --------------------------------------------------------------------------- */

template<>
std::vector<std::vector<PetscInt>> AF<honeycomb2D>::construct_sublat() {

  std::vector<std::vector<PetscInt>> sublat;
  short int tagx, tagy;
  PetscInt current;
  
  for (PetscInt s = 0; s < size; ++s)
    sublat.push_back(std::vector<PetscInt>{});

  tagx = 0;
  for (PetscInt ix = 0; ix < this->lattice.nspins_x; ++ix) {
    tagy = tagx;
    for (PetscInt iy = 0; iy < this->lattice.nspins_y; ++iy) {
      current = ix * this->lattice.nspins_y + iy;
      sublat[(tagy%2) ? 1 : 0].push_back(current);
      tagy ^= 1;
    }
    tagx ^= 1;
  }
  return sublat;
}

template<>
AF<honeycomb2D>::AF(honeycomb2D& m_lattice,
		    PetscInt m_size)
  : Sublattice<honeycomb2D>{m_lattice},
  nspins{m_lattice.num_spins()},
  size{2}
  {sublat = construct_sublat();}

/* --------------------------------------------------------------------------- */

template<>
std::vector<std::vector<PetscInt>> Striped<honeycomb2D>::construct_sublat() {

  std::vector<std::vector<PetscInt>> sublat;
  PetscInt current;
  short int tagx, tagy;
  
  for (PetscInt s = 0; s < size; ++s)
    sublat.push_back(std::vector<PetscInt>{});

  tagx = 0;
  for (PetscInt ix = 0; ix < this->lattice.nspins_x; ++ix) {
    tagy = 0;
    for (PetscInt iy = 0; iy < this->lattice.nspins_y; ++iy) {
      current = ix * this->lattice.nspins_y + iy;
      sublat[(tagy%2) ? 1 : 0].push_back(current);
      tagy ^= 1;
    }
  }

  /* MISSING THE OTHER FOUR SUBLATTICES.
     THE DIMENSIONS SHOULD BE MULTIPLE OF 4, IF NOT
     PERIODIC BOUNDARY CONDITIONS ARE NOT SATISFIED.
   */
  
  return sublat;
}

template<>
Striped<honeycomb2D>::Striped(honeycomb2D& m_lattice,
			      PetscInt m_size)
  : Sublattice<honeycomb2D>{m_lattice},
  nspins{m_lattice.nspins},
  size{6}
{sublat = construct_sublat();}
