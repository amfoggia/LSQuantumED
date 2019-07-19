#include "meson_config.h"
#include "random_disorder.hpp"

/* --------------------------------------------------------------------------- */
/* ------------------------ RandomDisorder Class ----------------------------- */
/* --------------------------------------------------------------------------- */

Tools::RandomDisorder::RandomDisorder(PetscInt m_nspins,
				      PetscBool m_reprod,
				      PetscReal m_min,
				      PetscReal m_max,
				      PetscInt m_rep)
  : min{m_min},
    max{m_max},
    rep{m_rep},
    nspins{m_nspins}
{
  PetscCalloc1(rep, &srand_seq);
  PetscCalloc1(nspins, &hi);

  srand_seq = srand_seq_init(m_reprod);
}

/* --------------------------------------------------------------------------- */

Tools::RandomDisorder::~RandomDisorder() {
  PetscFree(srand_seq);
  PetscFree(hi);
}

/* --------------------------------------------------------------------------- */

PetscInt* Tools::RandomDisorder::srand_seq_init(PetscBool reprod) {

  std::uniform_int_distribution<PetscInt> dist (0,1000);
  std::random_device rd;
  if (!reprod) {
    std::mt19937_64 gen_rand(rd());
    for (PetscInt i = 0; i < rep; ++i)
      srand_seq[i] = dist(gen_rand);
  }
    else {
      std::mt19937_64 gen_reprod(1);
      for (PetscInt i = 0; i < rep; ++i)
	srand_seq[i] = dist(gen_reprod);
    }

  return srand_seq;
}

/* --------------------------------------------------------------------------- */

void Tools::RandomDisorder::hi_init(PetscInt dis_index) {

  std::mt19937_64 gen(srand_seq[dis_index]);
  std::uniform_real_distribution<PetscReal> dist (min, max);
  for (PetscInt spin = 0; spin < nspins; ++spin)
    hi[spin] = dist(gen);
}

/* --------------------------------------------------------------------------- */
