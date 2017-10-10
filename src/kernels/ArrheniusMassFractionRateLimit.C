/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ArrheniusMassFractionRateLimit.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<ArrheniusMassFractionRateLimit>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Chemical reaction mass fraction decrease rate: Arrhenius equation with rate limit");
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
  return params;
}

ArrheniusMassFractionRateLimit::ArrheniusMassFractionRateLimit(const InputParameters & parameters)
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
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_name))
{
}

Real
ArrheniusMassFractionRateLimit::computeQpResidual()
{
  return _test[_i][_qp] * _mass_fraction_rate[_qp];
}

Real
ArrheniusMassFractionRateLimit::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] * _dmass_fraction_rate_dmass_fraction[_qp];
}

Real
ArrheniusMassFractionRateLimit::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _temperature_var) {
      return _test[_i][_qp] * _phi[_j][_qp] * _dmass_fraction_rate_dtemperature[_qp];
  }
  return 0.0;
}
