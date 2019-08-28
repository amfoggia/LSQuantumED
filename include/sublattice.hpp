#ifndef __SUBLATTICE_H
#define __SUBLATTICE_H

#include "meson_config.h"
#include "lattice.hpp"
#include <petscsys.h>
#include <vector>
#include <iostream>
#include <memory>

/* --------------------------------------------------------------------------- */
// --------------------------- Sublattice class ------------------------------ //
/* --------------------------------------------------------------------------- */

/**
 * @class Sublattice
 * @brief Creates list of spins that correspond to a certain sublattice type.
 * For example, an antiferromagnetic sublattice or a striped sublattice.
 * @tparam L Lattice type.
 */

template<typename L>
class Sublattice {

protected:
  L& lattice; /**< Lattice object. */
  std::vector<std::vector<PetscInt>> sublat; /**< Vector with sublattice(s). */

public:

  /**
   * @fn Sublattice(L&)
   * @brief Constructor.
   * @param[in] m_lattice Lattice object.
   */
  Sublattice(L& m_lattice)
    : lattice{m_lattice} {}

  /**
   * @fn construct_sublat
   * @brief Constructs the sublattice(s).
   * This is the function the user should provide if it wants to add a new sublattice derived class.
   * @return A vector with the sublattice(s).
   */
  virtual std::vector<std::vector<PetscInt>> construct_sublat() = 0;

  /**
   * @fn get_size
   * @brief Returns the number of sublattices created.
   * @return Number of sublattices.
   */
  virtual PetscInt get_size() const = 0;

  /**
   * @fn get_sl
   * @brief Returns a reference to the sublattice(s) vector.
   * @param sublat_i The sublattice you want.
   * @return Sublattice.
   */
  virtual const std::vector<PetscInt>& get_sl(PetscInt sublat_i) = 0;

  /**
   * @fn ~Sublattice
   * @brief Virtual destructor.
   */
  virtual ~Sublattice() {}
};

/* --------------------------------------------------------------------------- */
// ---------------------------------- AF class ------------------------------- //
/* --------------------------------------------------------------------------- */

/**
 * @class AF
 * @brief Creates two sublattices that correspond to an antiferromagnetic ground state.
 * @tparam L Lattice type.
 */
template<typename L>
class AF : public Sublattice<L> {

private:
  std::vector<std::vector<PetscInt>> sublat; /**< Vector with sublattice(s). */
  PetscInt nspins; /**< Number of spins in the system. */
  PetscInt size; /**< Number of sublattices of this type. */
  
  std::vector<std::vector<PetscInt>> construct_sublat() override;
  
public:

  /**
   * @fn AF(L&)
   * @brief Constructor.
   * @param[in] m_lattice Lattice object.
   * @param[in] m_size Number of sublattices of this type.
   */
  AF(L& m_lattice,
     PetscInt m_size);

  PetscInt get_size() const override {return size;};

  const std::vector<PetscInt>& get_sl(PetscInt sublat_i) override {return sublat[sublat_i];}
};

/* --------------------------------------------------------------------------- */
// ------------------------------ Striped class ------------------------------ //
/* --------------------------------------------------------------------------- */

/**
 * @class Striped
 * @brief Creates two sublattices that correspond to a striped ground state.
 * @tparam size Number of sublattices of this kind.
 * @tparam L Lattice type.
 */
template<typename L>
class Striped : public Sublattice<L> {

private:
  std::vector<std::vector<PetscInt>> sublat; /**< Vector with sublattice(s). */
  PetscInt nspins; /**< Number of spins in the system. */
  PetscInt size; /**< Number of sublattices of this type. */
  
  std::vector<std::vector<PetscInt>> construct_sublat() override;
public:

  /**
   * @fn Striped(L&)
   * @brief Constructor.
   * @param[in] m_lattice Lattice object.
   * @param[in] m_size Number of sublattices of this type.
   */
  Striped(L& m_lattice,
	  PetscInt m_size);
  
  PetscInt get_size() const override {return size;};

  const std::vector<PetscInt>& get_sl(PetscInt sublat_i) override {return sublat[sublat_i];}
};

#endif
