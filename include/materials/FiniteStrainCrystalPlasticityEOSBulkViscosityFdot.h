/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYEOSBULKVISCOSITYFDOT_H
#define FINITESTRAINCRYSTALPLASTICITYEOSBULKVISCOSITYFDOT_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityEOSBulkViscosity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * Bulk viscosity for oscillation damping is introduced (TrD = TrFdot)
 * Birch-Murnaghan equation of state is introduced
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class FiniteStrainCrystalPlasticityEOSBulkViscosityFdot;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityEOSBulkViscosityFdot>();

class FiniteStrainCrystalPlasticityEOSBulkViscosityFdot : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityEOSBulkViscosityFdot(const InputParameters & parameters);

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

#endif //FINITESTRAINCRYSTALPLASTICITYEOSBULKVISCOSITYFDOT_H
