#include "basis.hpp"
#include "boost/math/special_functions/binomial.hpp"
#include <stdexcept>

/* --------------------------------------------------------------------------- */
// -------------------------- Helper functions ------------------------------- //
/* --------------------------------------------------------------------------- */

PetscInt BasisHelper::search_elem(PetscInt index,
				  PetscInt nspins,
				  PetscInt nspins_up) {

  PetscInt pi = nspins - 1; // Position of i-eth one in elem
  PetscInt rest = index;    // Rest of (pi over i)
  PetscInt elem = 0;        // Number we obtain from the index
  PetscInt sup = nspins_up; // Current number one being treated
  PetscReal binomial;
  
  rest = index;
  pi = nspins - 1;
  elem = 0;
  sup = nspins_up;
  
  while (sup > 0) {

    binomial = (pi < sup ? 0 : boost::math::binomial_coefficient<double>((double)pi,(double)sup) );

    if ((PetscInt)std::floor(binomial + 0.5) <= rest) {
      rest -= binomial;
      elem |= (1ull << pi);
      sup = sup - 1;
    }
    pi -= 1;
  }
  return elem;
}

PetscInt Basis::compute_size() {

  PetscReal size = 1.0;
  for(PetscInt i = 1; i <= (nspins - nspins_up); ++i)
    size *= (static_cast<PetscReal> (i + nspins_up) / static_cast<PetscReal> (i));
  
  return floor(size + 0.5);
  // return PetscInt{size + 0.5};
}

/* --------------------------------------------------------------------------- */
// --------------------------- Constructors ---------------------------------- //
/* --------------------------------------------------------------------------- */

Basis::Basis(){}

/* --------------------------------------------------------------------------- */

Basis::Basis(Environment& env, PetscInt m_total_mag)
  : nspins{env.nspins},
    mpi_size{env.mpi_size},
    mpi_rank{env.mpi_rank}
{

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "Basis"};
#endif
    if (m_total_mag > nspins/2 || m_total_mag < -nspins/2)
      throw std::invalid_argument("The required constant magnetization has to be <= nspins/2.");

    max_mag = nspins/2;
    total_mag = m_total_mag;
    nspins_up = total_mag + nspins/2;
    size = compute_size();

    // Special case: there's not enough elems to split in processes
    if (size < mpi_size)
      throw std::invalid_argument("The size of the basis has to be at least equal to the number of MPI processes used.");

    PetscInt rest = size%mpi_size;
    local_size = size/mpi_size + (rest > mpi_rank);
    global_start_index = mpi_rank*local_size + (rest <= mpi_rank)*rest;

    // Initialize basis
    int_basis = generate_int_basis();
    int_basis.shrink_to_fit();

#ifdef DEVEL
    bit_basis = generate_bit_basis();
    bit_basis.shrink_to_fit();
#endif
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("Basis");
#endif
}

/* --------------------------------------------------------------------------- */

Basis::Basis(const Basis& rhs)
  : max_mag{rhs.max_mag},
    nspins{rhs.nspins},
    total_mag{rhs.total_mag},
    nspins_up{rhs.nspins_up},
    size{rhs.size},
    mpi_size{rhs.mpi_size},
    mpi_rank{rhs.mpi_rank},
    local_size{rhs.local_size},
    global_start_index{rhs.global_start_index},
    int_basis{rhs.int_basis}
#ifdef DEVEL
  ,bit_basis{rhs.bit_basis}
#endif
{}

/* --------------------------------------------------------------------------- */

Basis::Basis(Basis&& rhs) noexcept
  : max_mag{rhs.max_mag},
    nspins{rhs.nspins},
    total_mag{rhs.total_mag},
    nspins_up{rhs.nspins_up},
    size{rhs.size},
    mpi_size{rhs.mpi_size},
    mpi_rank{rhs.mpi_rank},
    local_size{rhs.local_size},
    global_start_index{rhs.global_start_index},
    int_basis{std::move(rhs.int_basis)}
#ifdef DEVEL
  ,bit_basis{std::move(rhs.bit_basis)}
#endif
{}

/* --------------------------------------------------------------------------- */
// ---------------------------- Assignments ---------------------------------- //
/* --------------------------------------------------------------------------- */

Basis& Basis::operator=(const Basis& rhs) {

  if (&rhs == this)
    return *this;

  max_mag = rhs.max_mag;
  nspins = rhs.nspins;
  total_mag = rhs.total_mag;
  nspins_up = rhs.nspins_up;
  size = rhs.size;
  mpi_size = rhs.mpi_size;
  mpi_rank = rhs.mpi_rank;
  local_size = rhs.local_size;
  global_start_index = rhs.global_start_index;
  int_basis = rhs.int_basis;
  int_basis.shrink_to_fit();

#ifdef DEVEL
  bit_basis = rhs.bit_basis;
  bit_basis.shrink_to_fit();
#endif
  
  return *this;
}

/* --------------------------------------------------------------------------- */

Basis& Basis::operator=(Basis&& rhs) noexcept {

  if (&rhs == this)
    return *this;

  max_mag = rhs.max_mag;
  nspins = rhs.nspins;
  total_mag = rhs.total_mag;
  nspins_up = rhs.nspins_up;
  size = rhs.size;
  mpi_size = rhs.mpi_size;
  mpi_rank = rhs.mpi_rank;
  local_size = rhs.local_size;
  global_start_index = rhs.global_start_index;
  int_basis = std::move(rhs.int_basis);
  int_basis.shrink_to_fit();

#ifdef DEVEL
  bit_basis = std::move(rhs.bit_basis);
  bit_basis.shrink_to_fit();
#endif
  
  return *this;
}

/* --------------------------------------------------------------------------- */
// ------------------------------ Generate basis ----------------------------- //
/* --------------------------------------------------------------------------- */

std::vector<PetscInt>& Basis::generate_int_basis() {
  // MAYBE THIS RETURNS A CONST OUTPUT

  // Check for the smallest int you can generate
  // with nspins_up in nspins
  PetscInt smallest_int = 0, local_smallest_int;
  for (PetscInt i = 0; i < nspins_up; ++i)
    smallest_int += 1<<i;

  // Get the first element of int_basis in
  // each process: local_smallest_int
  local_smallest_int = BasisHelper::search_elem(global_start_index, nspins, nspins_up);
  
  // Generate the basis with Gosper's hack
  PetscInt elem = local_smallest_int;
  int_basis.push_back(elem);
  for (PetscInt i = 1; i < local_size; ++i) {
    PetscInt u = elem & -elem;
    PetscInt v = u + elem;
    elem = (((v^elem)/u)>>2) + v;
    int_basis.push_back(elem);
  }
  return int_basis;
}

/* --------------------------------------------------------------------------- */

#ifdef DEVEL
std::vector<boost::dynamic_bitset<>>& Basis::generate_bit_basis() {

  for (PetscInt i = 0; i < local_size; ++i)
    bit_basis.push_back(boost::dynamic_bitset<>(nspins, int_basis[i]));
  return bit_basis;
}
#endif
