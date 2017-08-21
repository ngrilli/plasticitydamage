/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ArrheniusHeatEnergy.h"
#include "MooseMesh.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<ArrheniusHeatEnergy>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Chemical reaction heat energy density: Arrhenius equation");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addRequiredCoupledVar("mass_fraction",
                               "Mass fraction of the reactant");
  params.addParam<Real>("exponential_prefactor", 0.0, "Exponential prefactor in Arrhenius");
  params.addParam<Real>("exponential_coefficient", 0.0, "Exponential coefficient in Arrhenius");
  params.addParam<Real>("exponential_factor", 0.0, "Exponential factor in Arrhenius");
  return params;
}

ArrheniusHeatEnergy::ArrheniusHeatEnergy(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _mass_fraction(coupledValue("mass_fraction")),
    _mass_fraction_var(coupled("mass_fraction")),
    _exponential_prefactor(getParam<Real>("exponential_prefactor")),
    _exponential_coefficient(getParam<Real>("exponential_coefficient")),
    _exponential_factor(getParam<Real>("exponential_factor"))
{
}

Real
ArrheniusHeatEnergy::computeQpResidual()
{
  const Real & temperature = _u[_qp];

  return -_test[_i][_qp] * _mass_fraction[_qp] * _exponential_prefactor
         * std::exp(_exponential_coefficient - _exponential_factor / temperature);
}

Real
ArrheniusHeatEnergy::computeQpJacobian()
{
  const Real & temperature = _u[_qp];

  return -_test[_i][_qp] * _mass_fraction[_qp] * _exponential_prefactor
         * (_exponential_factor / (temperature * temperature)) * _phi[_j][_qp]
         * std::exp(_exponential_coefficient - _exponential_factor / temperature);
}

Real
ArrheniusHeatEnergy::computeQpOffDiagJacobian(unsigned int jvar)
{
  const Real & temperature = _u[_qp];

  if (jvar == _mass_fraction_var)
      return -_test[_i][_qp] * _phi[_j][_qp] * _exponential_prefactor
             * std::exp(_exponential_coefficient - _exponential_factor / temperature);

  return 0.0;
}
