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
 * @tparam size Number of sublattices of this kind.
 * @tparam L Lattice type.
 */

template<int size,typename L>
class Sublattice {

protected:
  L& lattice; /**< Lattice object. */

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
   * @return An array with the sublattice(s).
   */
  virtual std::array<std::vector<PetscInt>,size> construct_sublat() = 0;

  /**
   * @fn get_size
   * @brief Returns the number of sublattices created.
   * @return Number of sublattices.
   */
  virtual PetscInt get_size() const = 0;

  /**
   * @fn get_sl
   * @brief Returns a reference to the sublattice(s) array.
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
 * @tparam size Number of sublattices of this kind.
 * @tparam L Lattice type.
 */
template<int size,typename L>
class AF : public Sublattice<size,L> {

private:
  std::array<std::vector<PetscInt>,size> sublat; /**< Array with sublattice(s). */
  PetscInt nspins; /**< Number of spins in the system. */
  
  std::array<std::vector<PetscInt>,size> construct_sublat() override;
  
public:

  /**
   * @fn AF(L&)
   * @brief Constructor.
   * @param[in] m_lattice Lattice object.
   */
  AF(L& m_lattice)
    : Sublattice<size,L>{m_lattice},
    nspins{m_lattice.num_spins()}
  {sublat = construct_sublat();}

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
template<int size,typename L>
class Striped : public Sublattice<size,L> {

private:
  std::array<std::vector<PetscInt>,size> sublat; /**< Array with sublattice(s). */
  PetscInt nspins; /**< Number of spins in the system. */
  
  std::array<std::vector<PetscInt>,size> construct_sublat() override;
public:

  /**
   * @fn Striped(L&)
   * @brief Constructor.
   * @param[in] m_lattice Lattice object.
   */
  Striped(L& m_lattice)
    : Sublattice<size,L>{m_lattice},
    nspins{m_lattice.nspins}
  {sublat = construct_sublat();}

  PetscInt get_size() const override {return size;};

  const std::vector<PetscInt>& get_sl(PetscInt sublat_i) override {return sublat[sublat_i];}
};

#endif
