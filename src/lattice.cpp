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
    r = construct_r();
    r.shrink_to_fit();
    q = construct_q();
    q.shrink_to_fit();
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
  }
  return nnn;
}

/* --------------------------------------------------------------------------- */

std::vector<std::array<PetscInt,2>>& chain1D::construct_r() {

  for (PetscInt i = 0; i < nspins; ++i)
    r.push_back(std::array<PetscInt,2>{i,0});
  
  return r;
}

/* --------------------------------------------------------------------------- */

std::vector<std::array<PetscReal,2>>& chain1D::construct_q() {

  const PetscReal pi = boost::math::constants::pi<PetscReal>();
  PetscReal factor = 2.0 * pi / PetscReal(nspins);
  for (PetscInt i = 0; i < nspins; ++i)
    q.push_back(std::array<PetscReal,2>{PetscReal(i)*factor,0});

  return q;
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
    throw std::invalid_argument("The product between the number of spins in x direction \
                                   and the numer of spins in the y direction has to be equal \
                                   to the total number of spins.\n");

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "Lattice"};
#endif
    nn = get_nn();
    nn.shrink_to_fit();
    nnn = get_nnn();
    nnn.shrink_to_fit();
    r = construct_r();
    r.shrink_to_fit();
    q = construct_q();
    q.shrink_to_fit();
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
    }
  return nnn;
}

/* --------------------------------------------------------------------------- */

std::vector<std::array<PetscInt,2>>& square2D::construct_r() {
  
  for (PetscInt ix = 0; ix < nspins_x; ++ix)
    for (PetscInt iy = 0; iy < nspins_y; ++iy)
      r.push_back(std::array<PetscInt,2>{ix,iy});
  return r;
}

/* --------------------------------------------------------------------------- */

std::vector<std::array<PetscReal,2>>& square2D::construct_q() {

  const PetscReal pi = boost::math::constants::pi<PetscReal>();
  PetscReal factor_x = 2.0 * pi / PetscReal(nspins_x);
  PetscReal factor_y = 2.0 * pi / PetscReal(nspins_y);

  for (PetscInt ix = 0; ix < nspins_x; ++ix)
    for (PetscInt iy = 0; iy < nspins_y; ++iy)
      q.push_back(std::array<PetscReal,2>{PetscReal(ix)*factor_x, PetscReal(iy)*factor_y});
  
  return q;
}


