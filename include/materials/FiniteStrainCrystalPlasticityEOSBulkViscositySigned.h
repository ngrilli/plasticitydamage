/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYEOSBULKVISCOSITYSIGNED_H
#define FINITESTRAINCRYSTALPLASTICITYEOSBULKVISCOSITYSIGNED_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityEOSBulkViscosity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * Bulk viscosity for oscillation damping is introduced
 * Birch-Murnaghan equation of state is introduced
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class FiniteStrainCrystalPlasticityEOSBulkViscositySigned;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityEOSBulkViscositySigned>();

class FiniteStrainCrystalPlasticityEOSBulkViscositySigned : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityEOSBulkViscositySigned(const InputParameters & parameters);

protected:

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  // temperature
  const VariableValue & _temp;

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // Von Neuman coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  const Real _thermal_expansion;

  const Real _reference_temperature;

};

#endif //FINITESTRAINCRYSTALPLASTICITYEOSBULKVISCOSITYSIGNED_H
