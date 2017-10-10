/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeArrheniusMassFractionRateLimit.h"

template <>
InputParameters
validParams<ComputeArrheniusMassFractionRateLimit>()
{
  InputParameters params = validParams<Material>();
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addClassDescription("Arrhenius mass fraction rate from chemistry");
  params.addRequiredCoupledVar("mass_fraction","Mass fraction");
  params.addRequiredCoupledVar("temperature","Temperature");
  params.addParam<Real>("exponential_prefactor", 0.0, "Exponential prefactor in Arrhenius");
  params.addParam<Real>("exponential_coefficient", 0.0, "Exponential coefficient in Arrhenius");
  params.addParam<Real>("exponential_factor", 0.0, "Exponential factor in Arrhenius");
  params.addParam<Real>("rate_limit", 1.0e9, "Upper threshold for the reaction rate");
  return params;
}

ComputeArrheniusMassFractionRateLimit::ComputeArrheniusMassFractionRateLimit(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction(coupledValue("mass_fraction")),
    _mass_fraction_var(coupled("mass_fraction")),
    _mass_fraction_name(getVar("mass_fraction", 0)->name()),
    _temperature(coupledValue("temperature")),
    _temperature_var(coupled("temperature")),
    _temperature_name(getVar("temperature", 0)->name()),
    _exponential_prefactor(getParam<Real>("exponential_prefactor")),
    _exponential_coefficient(getParam<Real>("exponential_coefficient")),
    _exponential_factor(getParam<Real>("exponential_factor")),
    _rate_limit(getParam<Real>("rate_limit")),
    _mass_fraction_rate(declareProperty<Real>(_base_name + "mass_fraction_rate")),
    _dmass_fraction_rate_dtemperature(
        declarePropertyDerivative<Real>(_base_name + "mass_fraction_rate", _temperature_name)),
    _dmass_fraction_rate_dmass_fraction(
        declarePropertyDerivative<Real>(_base_name + "mass_fraction_rate", _mass_fraction_name))
{
}

void
ComputeArrheniusMassFractionRateLimit::computeQpProperties()
{
  Real mass_fraction = _mass_fraction[_qp];
  Real exp_reaction_rate;

  if (mass_fraction < 0.0) mass_fraction = 0.0; // check positive mass fraction

  exp_reaction_rate = _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);

  exp_reaction_rate = std::min(exp_reaction_rate,_rate_limit);

  _mass_fraction_rate[_qp] = mass_fraction * exp_reaction_rate;

  _dmass_fraction_rate_dmass_fraction[_qp] = 0.0;
  _dmass_fraction_rate_dtemperature[_qp] = 0.0;

  if (_fe_problem.currentlyComputingJacobian())
  {
    if (exp_reaction_rate < _rate_limit && mass_fraction > 0.0) {
      _dmass_fraction_rate_dmass_fraction[_qp] = exp_reaction_rate;
      _dmass_fraction_rate_dtemperature[_qp] = _mass_fraction_rate[_qp]
                                   * (_exponential_factor / (_temperature[_qp] * _temperature[_qp]));
    }
  }
}
