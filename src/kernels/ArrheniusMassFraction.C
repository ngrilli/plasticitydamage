/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ArrheniusMassFraction.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<ArrheniusMassFraction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Chemical reaction mass fraction decrease rate: Arrhenius equation");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("temperature",
                               "Temperature");
  params.addParam<Real>("exponential_prefactor", 0.0, "Exponential prefactor in Arrhenius");
  params.addParam<Real>("exponential_coefficient", 0.0, "Exponential coefficient in Arrhenius");
  params.addParam<Real>("exponential_factor", 0.0, "Exponential factor in Arrhenius");
  return params;
}

ArrheniusMassFraction::ArrheniusMassFraction(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _temperature(coupledValue("temperature")),
    _temperature_var(coupled("temperature")),
    _exponential_prefactor(getParam<Real>("exponential_prefactor")),
    _exponential_coefficient(getParam<Real>("exponential_coefficient")),
    _exponential_factor(getParam<Real>("exponential_factor"))
{
}

Real
ArrheniusMassFraction::computeQpResidual()
{
  const Real & mass_fraction = _u[_qp];

  return _test[_i][_qp] * mass_fraction * _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);
}

Real
ArrheniusMassFraction::computeQpJacobian()
{
  const Real & mass_fraction = _u[_qp];

  return _test[_i][_qp] * _phi[_j][_qp] * _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);
}

Real
ArrheniusMassFraction::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real & mass_fraction = _u[_qp];

  if (jvar == _temperature_var)
      return _test[_i][_qp] * mass_fraction * _exponential_prefactor
             * (_exponential_factor / (_temperature[_qp] * _temperature[_qp])) * _phi[_j][_qp]
             * std::exp(_exponential_coefficient - _exponential_factor / _temperature[_qp]);

  return 0.0;
}
