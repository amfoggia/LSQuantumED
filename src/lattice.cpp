#include "lattice.hpp"

/* --------------------------------------------------------------------------- */
// ------------------------------- 1D Chain Lattice -------------------------- //
/* --------------------------------------------------------------------------- */

chain1D::chain1D(Environment& env)
  : Lattice{env},
    type{lattice_type::chain1D}
{
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "Lattice"};
#endif
    nn = get_nn();
    nn.shrink_to_fit();
    nnn = get_nnn();
    nnn.shrink_to_fit();
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("Lattice");
#endif
}

/* --------------------------------------------------------------------------- */

std::vector<PetscInt>& chain1D::get_nn() {

  PetscInt neigh;

  for (PetscInt i = 0; i < nspins; ++i) {
    nn.push_back(i);
    neigh = (i+1)%nspins;
    nn.push_back(neigh);
    nnXsite.push_back(1);
  }
  return nn;
}

/* --------------------------------------------------------------------------- */

std::vector<PetscInt>& chain1D::get_nnn() {

  PetscInt neigh;

  for (PetscInt i = 0; i < nspins; ++i) {
    nnn.push_back(i);
    neigh = (i+2)%nspins;
    nnn.push_back(neigh);
    nnnXsite.push_back(1);
  }
  return nnn;
}

/* --------------------------------------------------------------------------- */
// ------------------------------ 2D Square Lattice -------------------------- //
/* --------------------------------------------------------------------------- */

square2D::square2D(Environment& env,
		   PetscInt m_nspins_x,
		   PetscInt m_nspins_y)
  : Lattice{env},
    nspins_x{m_nspins_x},
    nspins_y{m_nspins_y},
    type{lattice_type::square2D}
{
  if (m_nspins_x*m_nspins_y != env.nspins)
    throw std::invalid_argument("The product between the number of spins in x direction and the numer of spins in the y direction has to be equal to the total number of spins.");

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "Lattice"};
#endif
    nn = get_nn();
    nn.shrink_to_fit();
    nnn = get_nnn();
    nnn.shrink_to_fit();
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("Lattice");
#endif
}

/* --------------------------------------------------------------------------- */

std::vector<PetscInt>& square2D::get_nn() {

  PetscInt neighA, neighB;
  PetscInt ix1,iy1;
  int current_spin;

  for (PetscInt ix = 0; ix < nspins_x; ++ix)
    for (PetscInt iy = 0; iy < nspins_y; ++iy) {
      current_spin = ix * nspins_y + iy;
      nn.push_back(current_spin);

      iy1 = (iy+1)%nspins_y;
      neighA = ix * nspins_y + iy1;

      ix1 = (ix+1)%nspins_x;
      neighB = ix1 * nspins_y + iy;

      nn.push_back(neighA);
      nn.push_back(neighB);
      nnXsite.push_back(2);
    }
  return nn;
}

/* --------------------------------------------------------------------------- */

std::vector<PetscInt>& square2D::get_nnn() {

  PetscInt neighA, neighB;
  PetscInt ix1,iy1,iym1;
  int current_spin;

  for (PetscInt ix = 0; ix < nspins_x; ++ix)
    for (PetscInt iy = 0; iy < nspins_y; ++iy) {
      current_spin = ix * nspins_y + iy;
      nnn.push_back(current_spin);

      iy1 = (iy+1)%nspins_y;
      ix1 = (ix+1)%nspins_x;
      neighA = ix1 * nspins_y + iy1;

      iym1 = (iy-1+nspins_y)%nspins_y;
      ix1 = (ix+1)%nspins_x;
      neighB = ix1 * nspins_y + iym1;

      nnn.push_back(neighA);
      nnn.push_back(neighB);
      nnnXsite.push_back(2);
    }
  return nnn;
}

/* --------------------------------------------------------------------------- */
// ---------------------------- 2D Honeycomb Lattice ------------------------- //
/* --------------------------------------------------------------------------- */

honeycomb2D::honeycomb2D(Environment& env,
			 PetscInt m_nspins_x,
			 PetscInt m_nspins_y)
  : Lattice{env},
    nspins_x{m_nspins_x},
    nspins_y{m_nspins_y},
    type{lattice_type::honeycomb2D}
{
  if (m_nspins_x*m_nspins_y != env.nspins)
    throw std::invalid_argument("The product between the number of spins in x direction and the numer of spins in the y direction has to be equal to the total number of spins.");

  if (m_nspins_x%2 != 0 || m_nspins_y%2 != 0)
    throw std::invalid_argument("The number of spins in each direction has to be even.");
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "Lattice"};
#endif
    nn = get_nn();
    nn.shrink_to_fit();
    nnn = get_nnn();
    nnn.shrink_to_fit();
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("Lattice");
#endif
}

/* --------------------------------------------------------------------------- */

std::vector<PetscInt>& honeycomb2D::get_nn() {

  PetscInt neighA, neighB;
  PetscInt ix1,iy1;
  int current_spin;

  for (PetscInt ix = 0; ix < nspins_x; ++ix)
    for (PetscInt iy = 0; iy < nspins_y; ++iy) {

      current_spin = ix * nspins_y + iy;
      nn.push_back(current_spin);
      
      // Even rows of the lattice (xi even)
      if (ix%2 == 0) {

	if (current_spin%2 == 0) {
	  iy1 = (iy+1)%nspins_y;
	  neighA = ix * nspins_y + iy1;

	  ix1 = (ix+1)%nspins_x;
	  neighB = ix1 * nspins_y + iy;

	  nn.push_back(neighA);
	  nn.push_back(neighB);
	  nnXsite.push_back(2);
	}
	else {
	  iy1 = (iy+1)%nspins_y;
	  neighA = ix * nspins_y + iy1;
	  nn.push_back(neighA);
	  nn.push_back(-1);
	  nnXsite.push_back(1);
	}
      }

      // Odd rows of the lattice (xi odd)
      else {

	if (current_spin%2 != 0) {
	  iy1 = (iy+1)%nspins_y;
	  neighA = ix * nspins_y + iy1;

	  ix1 = (ix+1)%nspins_x;
	  neighB = ix1 * nspins_y + iy;

	  nn.push_back(neighA);
	  nn.push_back(neighB);
	  nnXsite.push_back(2);
	}
	else {
	  iy1 = (iy+1)%nspins_y;
	  neighA = ix * nspins_y + iy1;
	  nn.push_back(neighA);
	  nn.push_back(-1);
	  nnXsite.push_back(1);
	}
      }
    }
  return nn;
}

/* --------------------------------------------------------------------------- */

std::vector<PetscInt>& honeycomb2D::get_nnn() {

  PetscInt neighA, neighB, neighC;
  PetscInt ix1,iy1,iy2,iym1;
  int current_spin;

  for (PetscInt ix = 0; ix < nspins_x; ++ix)
    for (PetscInt iy = 0; iy < nspins_y; ++iy) {
      current_spin = ix * nspins_y + iy;
      nnn.push_back(current_spin);
      nnnXsite.push_back(3);

      iy2 = (iy+2)%nspins_y;
      neighA = ix * nspins_y + iy2;
      
      ix1 = (ix+1)%nspins_x;

      iy1 = (iy+1)%nspins_y;
      neighB = ix1 * nspins_y + iy1;

      iym1 = (iy-1+nspins_y)%nspins_y;
      neighC = ix1 * nspins_y + iym1;

      nnn.push_back(neighA);
      nnn.push_back(neighB);
      nnn.push_back(neighC);
    }
  
  return nnn;
}
