/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ArrheniusHeatEnergyRateLimit.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<ArrheniusHeatEnergyRateLimit>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Chemical reaction heat energy density: Arrhenius equation. Upper threshold for the reaction rate");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("mass_fraction",
                               "Mass fraction of the reactant");
  params.addRequiredCoupledVar("temperature",
                               "Temperature");
  params.addRequiredParam<MaterialPropertyName>(
      "mass_fraction_rate", "Decrease rate of the mass fraction due to chemistry");
  params.addParam<MaterialPropertyName>(
      "dmass_fraction_rate_dtemperature", "Material property name with derivative of mass_fraction_rate with temperature");
  params.addParam<MaterialPropertyName>(
      "dmass_fraction_rate_dmass_fraction", "Material property name with derivative of mass_fraction_rate with mass fraction");
  params.addParam<Real>("combustion_heat", 0.0, "Combustion heat");
  params.addRequiredParam<MaterialPropertyName>(
      "density", "Density");
  return params;
}

ArrheniusHeatEnergyRateLimit::ArrheniusHeatEnergyRateLimit(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction(coupledValue("mass_fraction")),
    _mass_fraction_var(coupled("mass_fraction")),
    _mass_fraction_name(getVar("mass_fraction", 0)->name()),
    _temperature(coupledValue("temperature")),
    _temperature_var(coupled("temperature")),
    _temperature_name(getVar("temperature", 0)->name()),
    _mass_fraction_rate(getMaterialProperty<Real>(_base_name + "mass_fraction_rate")),
    _dmass_fraction_rate_dtemperature(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _temperature_name)),
    _dmass_fraction_rate_dmass_fraction(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_name)),
    _combustion_heat(getParam<Real>("combustion_heat")),
    _density(getMaterialProperty<Real>(_base_name + "density"))
{
}

Real
ArrheniusHeatEnergyRateLimit::computeQpResidual()
{
  return -_test[_i][_qp] * _density[_qp] * _combustion_heat * _mass_fraction_rate[_qp];
}

Real
ArrheniusHeatEnergyRateLimit::computeQpJacobian()
{
  return -_test[_i][_qp] * _phi[_j][_qp] * _density[_qp] * _combustion_heat * _dmass_fraction_rate_dtemperature[_qp];
}

Real
ArrheniusHeatEnergyRateLimit::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _mass_fraction_var) {
    return -_test[_i][_qp] * _phi[_j][_qp] * _density[_qp] * _combustion_heat * _dmass_fraction_rate_dmass_fraction[_qp];
  }
  return 0.0;
}
