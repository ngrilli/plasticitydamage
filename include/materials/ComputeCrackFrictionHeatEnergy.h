/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTECRACKFRICTIONHEATENERGY_H
#define COMPUTECRACKFRICTIONHEATENERGY_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

/**
 * ComputeCrackFrictionHeatEnergy computes the energy from crack friction
 * and, if currentlyComputingJacobian, then the derivative of this quantity wrt total strain
 */
class ComputeCrackFrictionHeatEnergy : public DerivativeMaterialInterface<Material>
{
public:
  ComputeCrackFrictionHeatEnergy(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// optional parameter that allows multiple mechanics materials to be defined
  std::string _base_name;

  const VariableValue & _c;
  const VariableValue & _dcdx;
  const VariableValue & _dcdy;

  MaterialProperty<std::vector<Real>> & _crack_normal;
  MaterialProperty<Real> & _crack_normal_norm;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _strain_rate;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  const Real _friction_coefficient;
  const Real _l;

  MaterialProperty<std::vector<Real>> & _friction_force;
  MaterialProperty<std::vector<Real>> & _slide_velocity;

  MaterialProperty<Real> & _crack_surface_density;
  MaterialProperty<Real> & _heat_source_rate;

};

#endif // COMPUTECRACKFRICTIONHEATENERGY_H
