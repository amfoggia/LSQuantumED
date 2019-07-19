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
template<int size,typename L>
class AF;
template<int size,typename L>
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

  template<typename L>
  void get_elems(neigh_type nt,
		 Basis& basis,
		 PetscInt elem,
		 L& lattice,
		 PetscInt& ncol,
		 PetscInt& diag_elem,
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
  
  // /**
  //  * @fn construct_r
  //  * @brief Generates the list of distances of each spin to the origin.
  //  * @return Vector with list of distances of each spin to the origin.
  //  */
  // virtual std::vector<std::array<PetscInt,2>>& construct_r() = 0;

  // /**
  //  * @fn construct_q
  //  * @brief Generates the list of lattice vectors.
  //  * @return Vector with list of lattice vectors.
  //  */
  // virtual std::vector<std::array<PetscReal,2>>& construct_q() = 0;

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
   * @fn num_neighbours
   * @return Number of "effective" neighbours per lattice site.
   */
  virtual PetscInt num_neighbours() = 0;
  
  /**
   * @fn num_spins
   * @return Number of spins in the system
   */
  virtual PetscInt num_spins() {return nspins;};

  // /**
  //  * @fn get_r()
  //  * @brief Gives "read-only" access to the distance vector.
  //  * @return Distance vector.
  //  */
  // virtual const std::array<PetscInt,2>& get_r(PetscInt i) const = 0;
  
  // /**
  //  * @fn get_q()
  //  * @brief Gives "read-only" access to the lattice vectors.
  //  * @return Lattice vectors.
  //  */
  // virtual const std::array<PetscReal,2>& get_q(PetscInt i) const = 0;
  
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

  // std::vector<std::array<PetscInt,2>>& construct_r();
  // std::vector<std::array<PetscReal,2>>& construct_q();
  
  template<typename L> friend class Hamiltonian;
  template<int size,typename L> friend class AF;
  template<int size,typename L> friend class Striped;
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
  template<typename L> friend void HamiltHelper::get_elems(neigh_type nt,
							   Basis& basis,
							   PetscInt elem,
							   L& lattice,
							   PetscInt& ncol,
							   PetscInt& diag_elem,
							   PetscInt* coup_elems);

#ifdef DEVEL
public:
#endif
  std::vector<PetscInt> nn; /**< List of nearest neighbours of each spin. */  
  std::vector<PetscInt> nnn; /**< List of next-nearest neighbours of each spin. */
  // std::vector<std::array<PetscInt,2>> r; /**< List of distances of each spin to the origin. */
  // std::vector<std::array<PetscReal,2>> q; /**< List of lattice vectors. */

public:
  
  /**
   * @fn chain1D(Environment&)
   * @brief Constructor.
   * @param[in] env Environment object.
   */
  chain1D(Environment& env);

  lattice_type get_type() const {return type;}

  PetscInt num_neighbours() {return 1;}

  // const std::array<PetscInt,2>& get_r(PetscInt i) const {return r[i];}
  // const std::array<PetscReal,2>& get_q(PetscInt i) const {return q[i];}
};

/* --------------------------------------------------------------------------- */
// ------------------------------ 2D Square Lattice -------------------------- //
/* --------------------------------------------------------------------------- */

/**
 * @class square2D
 * @brief Creates a list of neighbours for the lattice sites for the 1D chain type of lattice.
 */

class square2D : public Lattice {

private:
  PetscInt nspins_x; /**< Number of spins in the x direction. */
  PetscInt nspins_y; /**< Number of spins in the y direction. */
  lattice_type type; /**< Type of lattice. */

  std::vector<PetscInt>& get_nn();

  std::vector<PetscInt>& get_nnn();

  // std::vector<std::array<PetscInt,2>>& construct_r();
  // std::vector<std::array<PetscReal,2>>& construct_q();
  
  template<typename L> friend class Hamiltonian;
  template<int size,typename L> friend class AF;
  template<int size,typename L> friend class Striped;
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
  template<typename L> friend void HamiltHelper::get_elems(neigh_type nt,
							   Basis& basis,
							   PetscInt elem,
							   L& lattice,
							   PetscInt& ncol,
							   PetscInt& diag_elem,
							   PetscInt* coup_elems);

#ifdef DEVEL
public:
#endif
  std::vector<PetscInt> nn; /**< List of nearest neighbours of each spin. */  
  std::vector<PetscInt> nnn; /**< List of next-nearest neighbours of each spin. */
  // std::vector<std::array<PetscInt,2>> r; /**< List of distances of each spin to the origin. */
  // std::vector<std::array<PetscReal,2>> q; /**< List of lattice vectors. */
  
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

  PetscInt num_neighbours() {return 2;}

  // const std::array<PetscInt,2>& get_r(PetscInt i) const {return r[i];}
  // const std::array<PetscReal,2>& get_q(PetscInt i) const {return q[i];}
};

#endif
