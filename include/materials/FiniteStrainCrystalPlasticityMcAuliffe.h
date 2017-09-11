/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYMCAULIFFE_H
#define FINITESTRAINCRYSTALPLASTICITYMCAULIFFE_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityMcAuliffe: McAuliffe and Waisman 2015
 */
class FiniteStrainCrystalPlasticityMcAuliffe;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityMcAuliffe>();

class FiniteStrainCrystalPlasticityMcAuliffe : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityMcAuliffe(const InputParameters & parameters);

protected:
  /**
   * This function set variables for internal variable solve.
   */
  virtual void preSolveStatevar();

  /**
   * This function solves internal variables.
   */
  virtual void solveStatevar();

  /**
   * This function update internal variable after solve.
   */
  virtual void postSolveStatevar();

  /**
   * Update elastic and plastic work
   */
  virtual void update_energies();

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  virtual void getSlipIncrements();

  const VariableValue & _c;

  /// Small number to avoid non-positive definiteness at or near complete damage
  Real _kdamage;

  MaterialProperty<Real> & _W0e;
  MaterialProperty<Real> & _W0p;
  MaterialProperty<Real> & _W0p_old;
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _dW0e_dstrain;
  MaterialProperty<RankTwoTensor> & _dW0p_dstrain;
  MaterialProperty<RankTwoTensor> & _pk2_undamaged;

  Real _W0p_tmp;
  Real _W0p_tmp_old;

};

#endif //FINITESTRAINCRYSTALPLASTICITYMCAULIFFE_H
