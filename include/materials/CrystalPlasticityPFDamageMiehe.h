/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CRYSTALPLASTICITYPFDAMAGEMIEHE_H
#define CRYSTALPLASTICITYPFDAMAGEMIEHE_H

#include "FiniteStrainCrystalPlasticity.h"

/**
 * FiniteStrainCrystalPlasticity uses the multiplicative decomposition of deformation gradient
 * and solves the PK2 stress residual equation at the intermediate configuration to evolve the material state.
 * The internal variables are updated using an interative predictor-corrector algorithm.
 * Backward Euler integration rule is used for the rate equations.
 */
class CrystalPlasticityPFDamageMiehe;

template<>
InputParameters validParams<CrystalPlasticityPFDamageMiehe>();

class CrystalPlasticityPFDamageMiehe : public FiniteStrainCrystalPlasticity
{
public:
  CrystalPlasticityPFDamageMiehe(const InputParameters & parameters);

protected:
  virtual void computeQpStress();
  virtual void preSolveQp();
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
  virtual void update_damage();

  /**
   * This function updates the slip increments.
   * And derivative of slip w.r.t. resolved shear stress.
   */
  virtual void getSlipIncrements();

  /**
   * This function calculate stress residual.
   */
  virtual void calcResidual( RankTwoTensor & );

  const MaterialProperty<Real> & _gc_prop;

  const VariableValue & _b;
  const VariableValue & _b_old;
  const VariableValue & _initialc;
  const VariableValue & _initialc_old;

  MaterialProperty<Real> & _c;
  MaterialProperty<Real> & _c_old;
  MaterialProperty<Real> & _damage;
  MaterialProperty<Real> & _damage_old;

  /// Small number to avoid non-positive definiteness at or near complete damage
  Real _kdamage;

  MaterialProperty<Real> & _x_bracket;
  MaterialProperty<Real> & _x_bracket_old;

  const Real _l;
  const Real _visco;
  const Real _Wc;

  // exponent in Murnaghan EOS
  const Real _n_Murnaghan;

  // reference bulk modulus
  const Real _Bulk_Modulus_Ref;

  // Von Neumann coefficient
  const Real _C0;

  // Landshoff coefficient
  const Real _C1;

  MaterialProperty<Real> & _W0e;
  MaterialProperty<Real> & _W0e_old;
  MaterialProperty<Real> & _W0p;
  MaterialProperty<Real> & _W0p_old;
  //MaterialProperty<RankTwoTensor> & _dstress_dc;
  MaterialProperty<RankTwoTensor> & _dW0e_dstrain;
  MaterialProperty<RankTwoTensor> & _dW0p_dstrain;
  MaterialProperty<RankTwoTensor> & _fe_old;
  MaterialProperty<RankTwoTensor> & _cauchy_out;

  Real _dfgrd_scale_factor_old;

  RankTwoTensor _dfgrd_tmp_old_substep;

  Real _W0e_tmp;
  Real _W0p_tmp;
  Real _damage_tmp;
  Real _b_tmp;
  Real _W0e_tmp_old;
  Real _W0p_tmp_old;
  Real _damage_tmp_old;
  Real _b_tmp_old;

  RankTwoTensor _fe_tmp_old;

};

#endif //CRYSTALPLASTICITYPFDAMAGEMIEHE_H
