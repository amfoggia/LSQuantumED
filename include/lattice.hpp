#ifndef __LATTICE_H
#define __LATTICE_H

#include "environment.hpp"
#include "meson_config.h"
#include <petscsys.h>
#include <vector>
#include <iostream>
#include <memory>
#include <boost/math/constants/constants.hpp>

/**
 * @enum lattice_type
 * @brief Type of lattice.
 */
enum class lattice_type {
  chain1D, /**< 1D chain type of lattice. */
  square2D, /**< 2D square type of lattice. */
  honeycomb2D /**< 2D honeycomb type of lattice. */
};

/**
 * @enum neigh_type
 * @brief Type of neighbour bond: nearest neighbours or next-nearest neighbours.
 */
enum class neigh_type {
  nn, /**< Nearest neighbours. */
  nnn, /**< Next-nearest neighbours. */
};

/* --------------------------------------------------------------------------- */
// ---------------------------------- Lattice -------------------------------- //
/* --------------------------------------------------------------------------- */

// Forward declatation of Basis class needed for the "friend" part
class Basis;

// Forward declatation of Hamiltonian class needed for the "friend" part
template<typename L>
class Hamiltonian;

// Forward declatation of AF and Striped class needed for the "friend" part
template<typename L>
class AF;
template<typename L>
class Striped;

// Forward declatation of HamiltHelper namespace needed for the "friend" part
namespace HamiltHelper {
  template<typename L>
  PetscInt get_elem_diag(neigh_type nt,
			 PetscInt basis_elem,
			 PetscInt nspins, 
			 L& lattice);

  template<typename L>
  void get_coup_elems(neigh_type nt,
		      Basis& basis,
		      PetscInt elem,
		      L& lattice,
		      PetscInt& ncol,
		      PetscInt* coup_elems);
}

/* --------------------------------------------------------------------------- */

/**
 * @class Lattice
 * @brief Creates a list of neighbours for the lattice sites depending on the type of lattice chosen.
 */

class Lattice {

protected:

  PetscInt nspins; /**< Number of spins in the system. */
  std::vector<PetscInt> nn; /**< List of nearest neighbours of each spin. */  
  std::vector<PetscInt> nnn; /**< List of next-nearest neighbours of each spin. */
  std::vector<PetscInt> nnXsite; /**< List with number of nearest neighbours per site. */
  std::vector<PetscInt> nnnXsite; /**< List with number of next-nearest neighbours per site. */
  
  /**
   * @fn get_nn
   * @brief Generates the list of next-nearest neighbours.
   * @return Vector with list of next-nearest neighbours for each lattice site.  
   */
  virtual std::vector<PetscInt>& get_nn() = 0;

  /**
   * @fn get_nnn
   * @brief Generates the list of next-nearest neighbours.
   * @return Vector with list of next-nearest neighbours for each lattice site.  
   */
  virtual std::vector<PetscInt>& get_nnn() = 0;

public:

  /**
   * @fn Lattice()
   * @brief Default constructor.
   */
  Lattice() = default;
  
  /**
   * @fn Lattice(Environment&)
   * @brief Constructor.
   * @param[in] env Environment object.
   */
  Lattice(Environment& env)
    : nspins{env.nspins} {}
  
  /**
   * @fn get_type
   * @return Lattice type.
   */
  virtual lattice_type get_type() const = 0;
  
  /**
   * @fn num_nn
   * @return Number of "effective" nearest neighbours per spin.
   */
  virtual PetscInt num_nn() = 0;

  /**
   * @fn num_nnn
   * @return Number of "effective" next-nearest neighbours per spin.
   */
  virtual PetscInt num_nnn() = 0;

  /**
   * @fn num_nnXsite(PetscInt)
   * @param[in] lat_site Lattice site.
   * @return Number of "effective" nearest neighbours of that lattice site.
   */
  virtual PetscInt num_nnXsite(PetscInt lat_site) = 0;

  /**
   * @fn num_nnnXsite(PetscInt)
   * @param[in] lat_site Lattice site.
   * @return Number of "effective" next-nearest neighbours of that lattice site.
   */
  virtual PetscInt num_nnnXsite(PetscInt lat_site) = 0;
  
  /**
   * @fn num_spins
   * @return Number of spins in the system
   */
  virtual PetscInt num_spins() {return nspins;};
  
  /**
   * @fn ~Lattice
   * @brief Virtual destructor.
   */
  virtual ~Lattice() {}
};

/* --------------------------------------------------------------------------- */
// ------------------------------- 1D Chain Lattice -------------------------- //
/* --------------------------------------------------------------------------- */

/**
 * @class chain1D
 * @brief Creates a list of neighbours for the lattice sites for the 1D chain type of lattice.
 */

class chain1D : public Lattice {

private:
  lattice_type type; /**< Type of lattice. */
  std::vector<PetscInt>& get_nn();
  std::vector<PetscInt>& get_nnn();

  template<typename L> friend class Hamiltonian;
  template<typename L> friend class AF;
  template<typename L> friend class Striped;
  template<typename L> friend PetscInt HamiltHelper::get_elem_diag(neigh_type nt,
								   PetscInt basis_elem,
								   PetscInt nspins, 
								   L& lattice);
  template<typename L> friend void HamiltHelper::get_coup_elems(neigh_type nt,
								Basis& basis,
								PetscInt elem,
								L& lattice,
								PetscInt& ncol,
								PetscInt* coup_elems);

#ifdef DEVEL
public:
#endif
  std::vector<PetscInt> nn; /**< List of nearest neighbours of each spin. */  
  std::vector<PetscInt> nnn; /**< List of next-nearest neighbours of each spin. */
  std::vector<PetscInt> nnXsite; /**< List with number of nearest neighbours per site. */
  std::vector<PetscInt> nnnXsite; /**< List with number of next-nearest neighbours per site. */

public:
  
  /**
   * @fn chain1D(Environment&)
   * @brief Constructor.
   * @param[in] env Environment object.
   */
  chain1D(Environment& env);

  lattice_type get_type() const {return type;}

  PetscInt num_nn() {return 1;}
  PetscInt num_nnn() {return 1;}
  PetscInt num_nnXsite(PetscInt lat_site) {return nnXsite[lat_site];}
  PetscInt num_nnnXsite(PetscInt lat_site) {return nnnXsite[lat_site];}
};

/* --------------------------------------------------------------------------- */
// ------------------------------ 2D Square Lattice -------------------------- //
/* --------------------------------------------------------------------------- */

/**
 * @class square2D
 * @brief Creates a list of neighbours for the lattice sites for the 2D square type of lattice.
 */

class square2D : public Lattice {

private:
  PetscInt nspins_x; /**< Number of spins in the x direction. */
  PetscInt nspins_y; /**< Number of spins in the y direction. */
  lattice_type type; /**< Type of lattice. */
  std::vector<PetscInt>& get_nn();
  std::vector<PetscInt>& get_nnn();
  
  template<typename L> friend class Hamiltonian;
  template<typename L> friend class AF;
  template<typename L> friend class Striped;
  template<typename L> friend PetscInt HamiltHelper::get_elem_diag(neigh_type nt,
								   PetscInt basis_elem,
								   PetscInt nspins, 
								   L& lattice);
  template<typename L> friend void HamiltHelper::get_coup_elems(neigh_type nt,
								Basis& basis,
								PetscInt elem,
								L& lattice,
								PetscInt& ncol,
								PetscInt* coup_elems);

#ifdef DEVEL
public:
#endif
  std::vector<PetscInt> nn; /**< List of nearest neighbours of each spin. */  
  std::vector<PetscInt> nnn; /**< List of next-nearest neighbours of each spin. */
  std::vector<PetscInt> nnXsite; /**< List with number of nearest neighbours per site. */
  std::vector<PetscInt> nnnXsite; /**< List with number of next-nearest neighbours per site. */
  
public:

  /**
   * @fn square2D(Environment&, PetscInt, PetscInt)
   * @brief Constructor.
   * @param[in] env Environment object.
   * @param[in] m_nspins_x Number of spins in x direction.
   * @param[in] m_nspins_y Number of spins in y direction.
   */
  square2D(Environment& env,
	   PetscInt m_nspins_x,
	   PetscInt m_nspins_y);
  
  lattice_type get_type() const {return type;}

  PetscInt num_nn() {return 2;}
  PetscInt num_nnn() {return 2;}
  PetscInt num_nnXsite(PetscInt lat_site) {return nnXsite[lat_site];}
  PetscInt num_nnnXsite(PetscInt lat_site) {return nnnXsite[lat_site];}
};

/* --------------------------------------------------------------------------- */
// ---------------------------- 2D Honeycomb Lattice ------------------------- //
/* --------------------------------------------------------------------------- */

/**
 * @class honeycomb2D
 * @brief Creates a list of neighbours for the lattice sites for the 2D honeycomb type of lattice.
 */

class honeycomb2D : public Lattice {

private:
  PetscInt nspins_x; /**< Number of spins in the x direction. */
  PetscInt nspins_y; /**< Number of spins in the y direction. */
  lattice_type type; /**< Type of lattice. */
  std::vector<PetscInt>& get_nn();
  std::vector<PetscInt>& get_nnn();

  template<typename L> friend class Hamiltonian;
  template<typename L> friend class AF;
  template<typename L> friend class Striped;
  template<typename L> friend PetscInt HamiltHelper::get_elem_diag(neigh_type nt,
  								   PetscInt basis_elem,
  								   PetscInt nspins, 
  								   L& lattice);
  template<typename L> friend void HamiltHelper::get_coup_elems(neigh_type nt,
  								Basis& basis,
  								PetscInt elem,
  								L& lattice,
  								PetscInt& ncol,
  								PetscInt* coup_elems);
  
#ifdef DEVEL
public:
#endif
  std::vector<PetscInt> nn; /**< List of nearest neighbours of each spin. */  
  std::vector<PetscInt> nnn; /**< List of next-nearest neighbours of each spin. */
  std::vector<PetscInt> nnXsite; /**< List with number of nearest neighbours per site. */
  std::vector<PetscInt> nnnXsite; /**< List with number of next-nearest neighbours per site. */

public:

  /**
   * @fn honeycomb2D(Environment&, PetscInt, PetscInt)
   * @brief Constructor.
   * @param[in] env Environment object.
   * @param[in] m_nspins_x Number of spins in x direction.
   * @param[in] m_nspins_y Number of spins in y direction.
   */
  honeycomb2D(Environment& env, PetscInt m_nspins_x, PetscInt m_nspins_y);
  lattice_type get_type() const {return type;}
  
  PetscInt num_nn() {return 2;}
  PetscInt num_nnn() {return 3;}
  PetscInt num_nnXsite(PetscInt lat_site) {return nnXsite[lat_site];}
  PetscInt num_nnnXsite(PetscInt lat_site) {return nnnXsite[lat_site];}
};

#endif
