/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ArrheniusMassFractionRateLimit3species.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<ArrheniusMassFractionRateLimit3species>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Chemical reaction mass fraction decrease rate: Arrhenius equation with rate limit for three species");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("mass_fraction_1","Mass fraction of the first specie");
  params.addRequiredCoupledVar("mass_fraction_2","Mass fraction of the second specie");
  params.addRequiredCoupledVar("mass_fraction_3","Mass fraction of the third specie");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addRequiredParam<MaterialPropertyName>(
    "mass_fraction_rate", "Decrease rate of the mass fraction due to chemistry");
  params.addParam<MaterialPropertyName>(
    "dmass_fraction_rate_dtemperature", "Material property name with derivative of mass_fraction_rate with temperature");
  params.addParam<MaterialPropertyName>(
    "dmass_fraction_rate_dmass_fraction_1", "Material property name with derivative of mass_fraction_rate with mass fraction 1");
  params.addParam<MaterialPropertyName>(
    "dmass_fraction_rate_dmass_fraction_2", "Material property name with derivative of mass_fraction_rate with mass fraction 2");
  params.addParam<MaterialPropertyName>(
    "dmass_fraction_rate_dmass_fraction_3", "Material property name with derivative of mass_fraction_rate with mass fraction 3");
  return params;
}

ArrheniusMassFractionRateLimit3species::ArrheniusMassFractionRateLimit3species(const InputParameters & parameters)
  : DerivativeMaterialInterface<Kernel>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction_1(coupledValue("mass_fraction_1")), // kernel variable
    _mass_fraction_1_var(coupled("mass_fraction_1")),
    _mass_fraction_1_name(getVar("mass_fraction_1", 0)->name()),
    _mass_fraction_2(coupledValue("mass_fraction_2")),
    _mass_fraction_2_var(coupled("mass_fraction_2")),
    _mass_fraction_2_name(getVar("mass_fraction_2", 0)->name()),
    _mass_fraction_3(coupledValue("mass_fraction_3")),
    _mass_fraction_3_var(coupled("mass_fraction_3")),
    _mass_fraction_3_name(getVar("mass_fraction_3", 0)->name()),
    _temperature(coupledValue("temperature")),
    _temperature_var(coupled("temperature")),
    _temperature_name(getVar("temperature", 0)->name()),
    _mass_fraction_rate(getMaterialProperty<Real>(_base_name + "mass_fraction_rate")),
    _dmass_fraction_rate_dtemperature(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _temperature_name)),
    _dmass_fraction_rate_dmass_fraction_1(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_1_name)),
    _dmass_fraction_rate_dmass_fraction_2(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_2_name)),
    _dmass_fraction_rate_dmass_fraction_3(
        getMaterialPropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_3_name))
{
}

Real
ArrheniusMassFractionRateLimit3species::computeQpResidual()
{
  return _test[_i][_qp] * _mass_fraction_rate[_qp];
}

Real
ArrheniusMassFractionRateLimit3species::computeQpJacobian()
{
  return _test[_i][_qp] * _phi[_j][_qp] * _dmass_fraction_rate_dmass_fraction_1[_qp];
}

Real
ArrheniusMassFractionRateLimit3species::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _temperature_var) {
      return _test[_i][_qp] * _phi[_j][_qp] * _dmass_fraction_rate_dtemperature[_qp];
  } else if (jvar == _mass_fraction_2_var) {
      return _test[_i][_qp] * _phi[_j][_qp] * _dmass_fraction_rate_dmass_fraction_2[_qp];
  } else if (jvar == _mass_fraction_3_var) {
      return _test[_i][_qp] * _phi[_j][_qp] * _dmass_fraction_rate_dmass_fraction_3[_qp];
  }
  return 0.0;
}
