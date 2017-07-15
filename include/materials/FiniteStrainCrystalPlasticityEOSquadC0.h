/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYEOSQUADC0_H
#define FINITESTRAINCRYSTALPLASTICITYEOSQUADC0_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityEOSquadC0 uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * Bulk viscosity for oscillation damping is introduced (quadratic C0)
 * Birch-Murnaghan equation of state is introduced
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class FiniteStrainCrystalPlasticityEOSquadC0;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityEOSquadC0>();

class FiniteStrainCrystalPlasticityEOSquadC0 : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityEOSquadC0(const InputParameters & parameters);

protected:

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // Von Neuman coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

};

#endif //FINITESTRAINCRYSTALPLASTICITYEOSQUADC0_H
