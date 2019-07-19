#ifndef __PHYS_H
#define __PHYS_H

#include "environment.hpp"
#include "hamiltonian.hpp"
#include "correlation.hpp"
#include "spinOPft.hpp"
#include "solver.hpp"
#include <sstream>


/**
 * @namespace Phys
 * @brief Functions to compute physical quantities.
 */

namespace Phys {

  /**
   * @struct DSF_data
   * @brief Necessary data for the computation of the dynamical structure factor.
   * @tparam Lattice type.
   */
  template<typename L>
  struct DSF_data {
    DSF_data() {}
    Basis* b0;                 /**< Basis of the Mag = 0 subspace. */
    Basis* b1;                 /**< Basis of the new subspace. */
    L* lat;                    /**< Lattice object. */
    Hamiltonian<L>* h0;        /**< Hamiltonian of the Mag = 0 subspace. */
    Hamiltonian<L>* h1;        /**< Hamiltonian of the new subspace. */
    PetscInt nev;              /**< Number of eigenvalues to compute. */
    PetscInt ncv;              /**< Maximum dimension of the subspace to be used by the solver. */
    PetscInt mpd;              /**< Maximum dimension allowed for the projected problem. */
    PetscBool disorder;        /**< Is there disorder in this computation? */
    PetscReal * hi;            /**< Disorder values. Pass NULL is there is no disorder realization. */
    const char * path;         /**< Path for the files regarding the computation of the DSF. */
  };

  /**
   * @fn SquaredSubLattMagAF(Environment&, Basis&, const std::array<std::vector<PetscInt>,2>&, Vec&)
   * @brief Computes the magnetic order parameter associated with the antiferromagnetic state.
   * Computes the antiferromagnetic squared sublattice magnetization.
   * @param[in] env Environment object.
   * @param[in] basis Basis object.
   * @param[in] sublat Array with sublattices configurations.
   * @param[in] state System state to compute the correlations.
   * @return Value of the parameter.
   */
  PetscScalar SquaredSubLattMagAF(Environment& env,
				  Basis& basis,
				  const std::array<std::vector<PetscInt>,2>& sublat,
				  Vec& state);

  /**
   * @fn SquaredSubLattMagSTR(Environment&, Basis&, const std::array<std::array<std::vector<PetscInt>,2>,2>&, Vec&)
   * @brief Computes the magnetic order parameter associated with the striped state.
   * Computes the striped squared sublattice magnetization.
   * @param[in] env Environment object.
   * @param[in] basis Basis object.
   * @param[in] sublat Array with sublattices configurations.
   * @param[in] state System state to compute the correlations.
   * @return Value of the parameter.
   */
  PetscScalar SquaredSubLattMagSTR(Environment& env,
				   Basis& basis,
				   const std::array<std::array<std::vector<PetscInt>,2>,2>& sublat,
				   Vec& state);

  /**
   * @fn DynStructFactor(Environment&, DSF_data<L>&, PetscScalar, Vec&, PetscInt)
   * @brief Computes the dynamical structure factor.
   * Computes the antiferromagnetic squared sublattice magnetization.
   * @tparam Lattice type.
   * @tparam Sq operator type.
   * @param[in] env Environment object.
   * @param[in] data DSF_data object with the data needed for this calculation.
   * @param[in] gs Ground state energy.
   * @param[in] state System state to compute the correlations.
   * @param[in] dis_iter Index of disorder realization iterations.
   * @return Value of the parameter.
   */
  template<typename L, template<typename> class SQ>
  PetscScalar DynStructFactor(Environment& env,
			      DSF_data<L>& data,
			      PetscScalar gs,
			      Vec& state,
			      PetscInt dis_iter = 0);
}

#include "phys.tmpl_func.hpp"

#endif
