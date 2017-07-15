/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYPHONON_H
#define FINITESTRAINCRYSTALPLASTICITYPHONON_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class FiniteStrainCrystalPlasticityPhonon;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityPhonon>();

class FiniteStrainCrystalPlasticityPhonon : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityPhonon(const InputParameters & parameters);

protected:
  /**
   * This function updates the slip increments.
   * And derivative of slip w.r.t. resolved shear stress.
   */
  virtual void getSlipIncrements();

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

#endif //FINITESTRAINCRYSTALPLASTICITYPHONON_H
