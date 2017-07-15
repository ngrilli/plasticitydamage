/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYBULKVISCOSITY_H
#define FINITESTRAINCRYSTALPLASTICITYBULKVISCOSITY_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityBulkViscosity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * Bulk viscosity for oscillation damping is introduced
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class FiniteStrainCrystalPlasticityBulkViscosity;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityBulkViscosity>();

class FiniteStrainCrystalPlasticityBulkViscosity : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityBulkViscosity(const InputParameters & parameters);

protected:

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  // Von Neuman coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

};

#endif //FINITESTRAINCRYSTALPLASTICITYBULKVISCOSITY_H
