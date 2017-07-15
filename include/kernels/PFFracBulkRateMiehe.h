/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFFRACBULKRATEMIEHE_H
#define PFFRACBULKRATEMIEHE_H

#include "KernelValue.h"
#include "RankTwoTensor.h"

// Forward Declarations
class PFFracBulkRateMiehe;

template <>
InputParameters validParams<PFFracBulkRateMiehe>();

/**
 * Phase field based fracture model
 * This kernel computes the residual and jacobian for bulk free energy contribution to c
 * Refer to Formulation: Miehe et. al., Int. J. Num. Methods Engg., 2010, 83. 1273-1311 Equation 63
 */
class PFFracBulkRateMiehe : public KernelValue
{
public:
  PFFracBulkRateMiehe(const InputParameters & parameters);

protected:
  virtual Real precomputeQpResidual();
  virtual Real precomputeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Critical energy release rate for fracture
  const MaterialProperty<Real> & _gc_prop;

  /// Contribution of umdamaged strain energy to damage evolution
  const MaterialProperty<Real> & _G0_pos;
  const MaterialProperty<Real> & _G0_pos_old;

  /// Variation of undamaged strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dG0_pos_dstrain;
  const MaterialProperty<RankTwoTensor> * _dG0_pos_dstrain_old;

  const MaterialProperty<Real> & _damage;
  const MaterialProperty<Real> & _damage_old;

  const MaterialProperty<Real> & _x_bracket;
  const MaterialProperty<Real> & _x_bracket_old;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  std::string _base_name;

  /// Characteristic length, controls damage zone thickness
  Real _l;

  /// Viscosity parameter ( visco -> 0, rate independent )
  Real _visco;
};

#endif // PFFRACBULKRATEMIEHE_H
