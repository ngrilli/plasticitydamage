/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYEOS_H
#define FINITESTRAINCRYSTALPLASTICITYEOS_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityEOS uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * Volumetric pressure from equation of state is introduced
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class FiniteStrainCrystalPlasticityEOS;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityEOS>();

class FiniteStrainCrystalPlasticityEOS : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityEOS(const InputParameters & parameters);

protected:

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

};

#endif //FINITESTRAINCRYSTALPLASTICITYEOS_H
