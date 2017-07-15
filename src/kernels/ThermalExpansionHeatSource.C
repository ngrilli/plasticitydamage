#include "ThermalExpansionHeatSource.h"

template <>
InputParameters
validParams<ThermalExpansionHeatSource>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Thermal expansion heat source kernel");
  params.addRequiredParam<Real>("thermal_expansion_heat","thermal expansion heat source coefficient: alpha E T_0 / (1-2 nu)");
  return params;
}

ThermalExpansionHeatSource::ThermalExpansionHeatSource(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _thermal_expansion_heat(getParam<Real>("thermal_expansion_heat")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
ThermalExpansionHeatSource::computeQpResidual()
{
  RankTwoTensor strain_rate;
  strain_rate = ( _mechanical_strain[_qp] - _mechanical_strain_old[_qp] ) / _dt;

  return _thermal_expansion_heat * strain_rate.trace() * _test[_i][_qp];
}

Real
ThermalExpansionHeatSource::computeQpJacobian()
{
  return -0.0 * _phi[_j][_qp] * _test[_i][_qp];
}

Real
ThermalExpansionHeatSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  unsigned int i, j, l, h;

  val = 0.0;
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k]) {
      val = _thermal_expansion_heat * _grad_phi[_j][_qp](k) * _test[_i][_qp] / _dt;
    }
  }
  return val;
}
