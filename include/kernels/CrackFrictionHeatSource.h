#ifndef CRACKFRICTIONHEATSOURCE_H
#define CRACKFRICTIONHEATSOURCE_H

#include "HeatSource.h"
#include "RankTwoTensor.h"

// Forward Declarations
class CrackFrictionHeatSource;

template <>
InputParameters validParams<CrackFrictionHeatSource>();

/**
 * This kernel calculates the heat source term corresponding to crack friction
 */
class CrackFrictionHeatSource : public HeatSource
{
public:
  CrackFrictionHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:

  std::string _base_name;

  const VariableValue & _dcdx;
  const VariableValue & _dcdy;
  //const VariableValue & _dcdz;

  const MaterialProperty<std::vector<Real>> & _crack_normal;
  const MaterialProperty<Real> & _crack_normal_norm;

  const MaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankTwoTensor> & _strain_rate;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  const Real _friction_coefficient;

  const MaterialProperty<std::vector<Real>> & _friction_force;
  const MaterialProperty<std::vector<Real>> & _slide_velocity;
  const MaterialProperty<Real> & _heat_source_rate;

  /// Coupled displacement variables
  const unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

};

#endif // CRACKFRICTIONHEATSOURCE_H
