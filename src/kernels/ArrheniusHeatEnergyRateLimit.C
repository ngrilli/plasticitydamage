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
  params.addParam<Real>("exponential_prefactor", 0.0, "Exponential prefactor in Arrhenius");
  params.addParam<Real>("exponential_coefficient", 0.0, "Exponential coefficient in Arrhenius");
  params.addParam<Real>("exponential_factor", 0.0, "Exponential factor in Arrhenius");
  params.addParam<Real>("rate_limit", 1.0e9, "Upper threshold for the reaction rate");
  return params;
}

ArrheniusHeatEnergyRateLimit::ArrheniusHeatEnergyRateLimit(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction(coupledValue("mass_fraction")),
    _mass_fraction_var(coupled("mass_fraction")),
    _exponential_prefactor(getParam<Real>("exponential_prefactor")),
    _exponential_coefficient(getParam<Real>("exponential_coefficient")),
    _exponential_factor(getParam<Real>("exponential_factor")),
    _rate_limit(getParam<Real>("rate_limit"))
{
}

Real
ArrheniusHeatEnergyRateLimit::computeQpResidual()
{
  const Real & temperature = _u[_qp];
  Real mass_fraction = _mass_fraction[_qp];
  Real reaction_rate;

  if (mass_fraction < 0.0) mass_fraction = 0.0; // check positive mass fraction

  reaction_rate = mass_fraction * _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / temperature);

  reaction_rate = std::max(reaction_rate,_rate_limit);

  return -_test[_i][_qp] * reaction_rate;
}

Real
ArrheniusHeatEnergyRateLimit::computeQpJacobian()
{
  const Real & temperature = _u[_qp];
  Real reaction_rate, exponential_rate;

  exponential_rate = _exponential_prefactor
                     * std::exp(_exponential_coefficient - _exponential_factor / temperature);

  reaction_rate = _mass_fraction[_qp] * exponential_rate;

  if (reaction_rate < _rate_limit && _mass_fraction[_qp] > 0.0) {
    return -_test[_i][_qp] * reaction_rate
         * (_exponential_factor / (temperature * temperature)) * _phi[_j][_qp];
  }
  else {
    return 0.0;
  }
}

Real
ArrheniusHeatEnergyRateLimit::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real & temperature = _u[_qp];
  Real reaction_rate, exponential_rate;

  exponential_rate = _exponential_prefactor
                     * std::exp(_exponential_coefficient - _exponential_factor / temperature);

  reaction_rate = _mass_fraction[_qp] * exponential_rate;

  if (jvar == _mass_fraction_var) {
    if (reaction_rate < _rate_limit && _mass_fraction[_qp] > 0.0) {
      return -_test[_i][_qp] * _phi[_j][_qp] * exponential_rate;
    }
    else {
      return 0.0;
    }
  }

  return 0.0;
}
