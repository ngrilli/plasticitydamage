/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef PFTHERMALCONDUCTIVITYCOMPRESSIVE_H
#define PFTHERMALCONDUCTIVITYCOMPRESSIVE_H

#include "Material.h"
#include "RankTwoTensor.h"

// Forward Declarations
class PFThermalConductivityCompressive;
class Function;

template <>
InputParameters validParams<PFThermalConductivityCompressive>();

/**
 * Effective thermal conductivity of damaged material
 */
class PFThermalConductivityCompressive : public Material
{
public:
  PFThermalConductivityCompressive(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const VariableValue & _c;

  const Real _k_m;

  const Real _k_c;

  MaterialProperty<Real> & _k;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<Real> & _friction_normal_force;
};

#endif // PFTHERMALCONDUCTIVITYCOMPRESSIVE_H
