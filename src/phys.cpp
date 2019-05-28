#include "phys.hpp"

/* --------------------------------------------------------------------------- */
// ------------------------ Magnetic Order Parameter ------------------------- //
/* --------------------------------------------------------------------------- */

PetscScalar Phys::SquaredSubLattMagAF(Environment& env,
				      Basis& basis,
				      const std::array<std::vector<PetscInt>,2>& sublat,
				      Vec& state) {

  PetscScalar mAF = 0.0;
  
#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "MagAF"};
#endif
    
    SiSj cij{env,basis};
    
    PetscScalar exp_val = 0.0;
    PetscInt spin_i, spin_j;
    
    for (PetscInt lat = 0; lat < 2; ++lat) {
      for (PetscInt i = 0; i < PetscInt(sublat[lat].size()); ++i) {
	spin_i = sublat[lat][i];
	for (PetscInt j = 0; j < PetscInt(sublat[lat].size()); ++j) {
	  spin_j = sublat[lat][j];
	  exp_val += cij.ExpectVal(env,spin_i, spin_j, state);
	}
      }
    }
    
    mAF = exp_val * 8.0 / PetscScalar(env.nspins) / (PetscScalar(env.nspins + 4.0));
    
#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("MagAF");
#endif
  
  return mAF;
}

/* --------------------------------------------------------------------------- */

PetscScalar Phys::SquaredSubLattMagSTR(Environment& env,
				       Basis& basis,
				       const std::array<std::array<std::vector<PetscInt>,2>,2>& sublat,
				       Vec& state) {

  PetscScalar mSTR = 0.0;

#ifdef TIME_CODE
  {
    Tools::ScopedTimer _timer_{env.tm, "MagSTR"};
#endif
    
    SiSj cij{env,basis};

    PetscScalar exp_val = 0.0;
    PetscInt spin_i, spin_j;

    for (PetscInt str = 0; str < 2; ++str)
      for (PetscInt lat = 0; lat < 2; ++lat)
	for (PetscInt i = 0; i < PetscInt(sublat[str][lat].size()); ++i) {
	  spin_i = sublat[str][lat][i];
	  for (PetscInt j = 0; j < PetscInt(sublat[str][lat].size()); ++j) {
	    spin_j = sublat[str][lat][j];
	    exp_val += cij.ExpectVal(env,spin_i, spin_j, state);
	  }
	}

    mSTR = exp_val * 4.0 / PetscScalar(env.nspins) / (PetscScalar(env.nspins + 4.0));

#ifdef TIME_CODE
  }
  env.tm.PrintTimeInfoFunc("MagSTR");
#endif
  
  return mSTR;
}
