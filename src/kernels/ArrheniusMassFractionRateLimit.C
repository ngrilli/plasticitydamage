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
  params.addRequiredCoupledVar("temperature",
                               "Temperature");
  params.addParam<Real>("exponential_prefactor", 0.0, "Exponential prefactor in Arrhenius");
  params.addParam<Real>("exponential_coefficient", 0.0, "Exponential coefficient in Arrhenius");
  params.addParam<Real>("exponential_factor", 0.0, "Exponential factor in Arrhenius");
  params.addParam<Real>("rate_limit", 1.0e9, "Upper threshold for the reaction rate");
  return params;
}

ArrheniusMassFractionRateLimit::ArrheniusMassFractionRateLimit(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _temperature(coupledValue("temperature")),
    _temperature_var(coupled("temperature")),
    _exponential_prefactor(getParam<Real>("exponential_prefactor")),
    _exponential_coefficient(getParam<Real>("exponential_coefficient")),
    _exponential_factor(getParam<Real>("exponential_factor")),
    _rate_limit(getParam<Real>("rate_limit"))
{
}

Real
ArrheniusMassFractionRateLimit::computeQpResidual()
{
  Real mass_fraction = _u[_qp];
  Real reaction_rate;

  if (mass_fraction < 0.0) mass_fraction = 0.0; // check positive mass fraction

  reaction_rate = mass_fraction * _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);

  reaction_rate = std::min(reaction_rate,_rate_limit);

  return _test[_i][_qp] * reaction_rate;
}

Real
ArrheniusMassFractionRateLimit::computeQpJacobian()
{
  const Real & mass_fraction = _u[_qp];
  Real reaction_rate, exponential_rate;

  exponential_rate = _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);

  reaction_rate = mass_fraction * exponential_rate;

  if (reaction_rate < _rate_limit && mass_fraction > 0.0) {
    return _test[_i][_qp] * _phi[_j][_qp] * exponential_rate;
  }
  else {
    return 0.0;
  }
}

Real
ArrheniusMassFractionRateLimit::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real & mass_fraction = _u[_qp];
  Real reaction_rate, exponential_rate;

  exponential_rate = _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);

  reaction_rate = mass_fraction * exponential_rate;

  if (jvar == _temperature_var) {
    if (reaction_rate < _rate_limit && mass_fraction > 0.0) {
      return _test[_i][_qp] * reaction_rate
             * (_exponential_factor / (_temperature[_qp] * _temperature[_qp])) * _phi[_j][_qp];
    }
    else {
      return 0.0;
    }
  }

  return 0.0;
}
