/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef FINITESTRAINCRYSTALPLASTICITYDAMAGEPRINCIPALSTRAINS_H
#define FINITESTRAINCRYSTALPLASTICITYDAMAGEPRINCIPALSTRAINS_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticityDamagePrincipalStrains uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 * Principal valeus of the Lagrangian strain used to calculate damage
 */
class FiniteStrainCrystalPlasticityDamagePrincipalStrains;

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityDamagePrincipalStrains>();

class FiniteStrainCrystalPlasticityDamagePrincipalStrains : public FiniteStrainCrystalPlasticity
{
public:
  FiniteStrainCrystalPlasticityDamagePrincipalStrains(const InputParameters & parameters);

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

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  MaterialProperty<Real> & _W0e;
  MaterialProperty<Real> & _W0p;
  MaterialProperty<Real> & _W0p_old;
  MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _dW0e_dstrain;
  MaterialProperty<RankTwoTensor> & _dW0p_dstrain;
  MaterialProperty<RankTwoTensor> & _pk2_undamaged;
  MaterialProperty<RankTwoTensor> & _fe_out; // Elastic deformation gradient for output

  Real _W0p_tmp;
  Real _W0p_tmp_old;

  std::vector<RankTwoTensor> _etens;
  std::vector<Real> _epos;
  std::vector<Real> _eigval;
  RankTwoTensor _eigvec;

};

#endif //FINITESTRAINCRYSTALPLASTICITYDAMAGEPRINCIPALSTRAINS_H
