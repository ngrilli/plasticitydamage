/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFFRACBULKRATENOBLIMITED_H
#define PFFRACBULKRATENOBLIMITED_H

#include "KernelValue.h"
#include "RankTwoTensor.h"

// Forward Declarations
class PFFracBulkRateNobLimited;

template <>
InputParameters validParams<PFFracBulkRateNobLimited>();

/**
 * Phase field based fracture model
 * This kernel computes the residual and jacobian for bulk free energy contribution to c
 * Refer to Formulation: Miehe et. al., Int. J. Num. Methods Engg., 2010, 83. 1273-1311 Equation 63
 */
class PFFracBulkRateNobLimited : public KernelValue
{
public:
  PFFracBulkRateNobLimited(const InputParameters & parameters);

protected:
  virtual Real precomputeQpResidual();
  virtual Real precomputeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Critical energy release rate for fracture
  const MaterialProperty<Real> & _gc_prop;

  /// Contribution of umdamaged strain energy to damage evolution
  const MaterialProperty<Real> & _G0_pos;

  /// Variation of undamaged strain energy driving damage evolution with strain
  const MaterialProperty<RankTwoTensor> * _dG0_pos_dstrain;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;
  std::string _base_name;
  const VariableValue & _d2c_dx2;
  const VariableValue & _d2c_dy2;
  const VariableValue & _d2c_dz2;

  /// Characteristic length, controls damage zone thickness
  Real _l;

  /// Viscosity parameter ( visco -> 0, rate independent )
  Real _visco;
};

#endif // PFFRACBULKRATENOBLIMITED_H
