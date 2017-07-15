/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef LINEARISOELASTICPFDAMAGEMIEHE_H
#define LINEARISOELASTICPFDAMAGEMIEHE_H

#include "ComputeStressBase.h"
#include "Function.h"

/**
 * Phase-field fracture
 * This class computes the energy contribution to damage growth
 * Small strain Isotropic Elastic formulation
 * Stiffness matrix scaled for heterogeneous elasticity property
 */
class LinearIsoElasticPFDamageMiehe : public ComputeStressBase
{
public:
  LinearIsoElasticPFDamageMiehe(const InputParameters & parameters);

protected:
  virtual void computeQpStress();
  virtual void updateVar();
  virtual void updateJacobian();

  /// Critical energy release rate for fracture
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

  MaterialProperty<Real> & _G0_pos;
  MaterialProperty<Real> & _G0_pos_old;
  MaterialProperty<Real> & _x_bracket;
  MaterialProperty<Real> & _x_bracket_old;

  /// Characteristic length, controls damage zone thickness
  Real _l;

  /// Viscosity parameter ( visco -> 0, rate independent )
  Real _visco;

  std::vector<RankTwoTensor> _etens;
  std::vector<Real> _epos;
  std::vector<Real> _eigval;
  RankTwoTensor _eigvec;

};

#endif // LINEARISOELASTICPFDAMAGEMIEHE_H
