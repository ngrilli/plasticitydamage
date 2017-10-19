/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "HeatConductionMaterial3species.h"
#include "Function.h"

// libmesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<HeatConductionMaterial3species>()
{
  InputParameters params = validParams<Material>();

  params.addCoupledVar("temp", "Coupled Temperature");
  params.addCoupledVar("mass_fraction_1", "Mass fraction of the first specie");
  params.addCoupledVar("mass_fraction_2", "Mass fraction of the second specie");
  params.addCoupledVar("mass_fraction_3", "Mass fraction of the third specie");

  params.addParam<Real>("thermal_conductivity_1", "The thermal conductivity value of the first specie");
  params.addParam<FunctionName>("thermal_conductivity_temperature_function_1",
                                "",
                                "Thermal conductivity as a function of temperature of the first specie.");
  params.addParam<Real>("thermal_conductivity_2", "The thermal conductivity value of the second specie");
  params.addParam<FunctionName>("thermal_conductivity_temperature_function_2",
                                "",
                                "Thermal conductivity as a function of temperature of the second specie.");
  params.addParam<Real>("thermal_conductivity_3", "The thermal conductivity value of the third specie");
  params.addParam<FunctionName>("thermal_conductivity_temperature_function_3",
                                "",
                                "Thermal conductivity as a function of temperature of the third specie.");

  params.addParam<Real>("specific_heat_1", "The specific heat value of the first specie");
  params.addParam<FunctionName>(
      "specific_heat_temperature_function_1", "", "Specific heat as a function of temperature of the first specie.");
  params.addParam<Real>("specific_heat_2", "The specific heat value of the second specie");
  params.addParam<FunctionName>(
      "specific_heat_temperature_function_2", "", "Specific heat as a function of temperature of the second specie.");
  params.addParam<Real>("specific_heat_3", "The specific heat value of the third specie");
  params.addParam<FunctionName>(
      "specific_heat_temperature_function_3", "", "Specific heat as a function of temperature of the third specie.");
  params.addClassDescription("Heat conduction for three species");

  return params;
}

HeatConductionMaterial3species::HeatConductionMaterial3species(const InputParameters & parameters)
  : Material(parameters),

    _has_temp(isCoupled("temp")),
    _temperature(_has_temp ? coupledValue("temp") : _zero),
    _has_mass_fraction_1(isCoupled("mass_fraction_1")),
    _mass_fraction_1(_has_mass_fraction_1 ? coupledValue("mass_fraction_1") : _zero),
    _has_mass_fraction_2(isCoupled("mass_fraction_2")),
    _mass_fraction_2(_has_mass_fraction_2 ? coupledValue("mass_fraction_2") : _zero),
    _has_mass_fraction_3(isCoupled("mass_fraction_3")),
    _mass_fraction_3(_has_mass_fraction_3 ? coupledValue("mass_fraction_3") : _zero),
    _my_thermal_conductivity_1(
        isParamValid("thermal_conductivity_1") ? getParam<Real>("thermal_conductivity_1") : 0),
    _my_specific_heat_1(isParamValid("specific_heat_1") ? getParam<Real>("specific_heat_1") : 0),
    _my_thermal_conductivity_2(
        isParamValid("thermal_conductivity_2") ? getParam<Real>("thermal_conductivity_2") : 0),
    _my_specific_heat_2(isParamValid("specific_heat_2") ? getParam<Real>("specific_heat_2") : 0),
    _my_thermal_conductivity_3(
        isParamValid("thermal_conductivity_3") ? getParam<Real>("thermal_conductivity_3") : 0),
    _my_specific_heat_3(isParamValid("specific_heat_3") ? getParam<Real>("specific_heat_3") : 0),

    _thermal_conductivity(declareProperty<Real>("thermal_conductivity")),
    _thermal_conductivity_dT(declareProperty<Real>("thermal_conductivity_dT")),
    _thermal_conductivity_temperature_function_1(
        getParam<FunctionName>("thermal_conductivity_temperature_function_1") != ""
            ? &getFunction("thermal_conductivity_temperature_function_1")
            : NULL),
    _thermal_conductivity_temperature_function_2(
        getParam<FunctionName>("thermal_conductivity_temperature_function_2") != ""
            ? &getFunction("thermal_conductivity_temperature_function_2")
            : NULL),
    _thermal_conductivity_temperature_function_3(
        getParam<FunctionName>("thermal_conductivity_temperature_function_3") != ""
            ? &getFunction("thermal_conductivity_temperature_function_3")
            : NULL),

    _specific_heat(declareProperty<Real>("specific_heat")),
    _specific_heat_temperature_function_1(
        getParam<FunctionName>("specific_heat_temperature_function_1") != ""
            ? &getFunction("specific_heat_temperature_function_1")
            : NULL),
    _specific_heat_temperature_function_2(
        getParam<FunctionName>("specific_heat_temperature_function_2") != ""
            ? &getFunction("specific_heat_temperature_function_2")
            : NULL),
    _specific_heat_temperature_function_3(
        getParam<FunctionName>("specific_heat_temperature_function_3") != ""
            ? &getFunction("specific_heat_temperature_function_3")
            : NULL)
{
  if (_thermal_conductivity_temperature_function_1 && !_has_temp)
  {
    mooseError("Must couple with temperature if using thermal conductivity function");
  }
  if (_thermal_conductivity_temperature_function_2 && !_has_temp)
  {
    mooseError("Must couple with temperature if using thermal conductivity function");
  }
  if (_thermal_conductivity_temperature_function_3 && !_has_temp)
  {
    mooseError("Must couple with temperature if using thermal conductivity function");
  }
  if (isParamValid("thermal_conductivity_1") && _thermal_conductivity_temperature_function_1)
  {
    mooseError(
        "Cannot define both thermal conductivity and thermal conductivity temperature function");
  }
  if (isParamValid("thermal_conductivity_2") && _thermal_conductivity_temperature_function_2)
  {
    mooseError(
        "Cannot define both thermal conductivity and thermal conductivity temperature function");
  }
  if (isParamValid("thermal_conductivity_3") && _thermal_conductivity_temperature_function_3)
  {
    mooseError(
        "Cannot define both thermal conductivity and thermal conductivity temperature function");
  }
  if (_specific_heat_temperature_function_1 && !_has_temp)
  {
    mooseError("Must couple with temperature if using specific heat function");
  }
  if (_specific_heat_temperature_function_2 && !_has_temp)
  {
    mooseError("Must couple with temperature if using specific heat function");
  }
  if (_specific_heat_temperature_function_3 && !_has_temp)
  {
    mooseError("Must couple with temperature if using specific heat function");
  }
  if (isParamValid("specific_heat_1") && _specific_heat_temperature_function_1)
  {
    mooseError("Cannot define both specific heat and specific heat temperature function");
  }
  if (isParamValid("specific_heat_2") && _specific_heat_temperature_function_2)
  {
    mooseError("Cannot define both specific heat and specific heat temperature function");
  }
  if (isParamValid("specific_heat_3") && _specific_heat_temperature_function_3)
  {
    mooseError("Cannot define both specific heat and specific heat temperature function");
  }
}

void
HeatConductionMaterial3species::computeProperties()
{
  for (unsigned int qp(0); qp < _qrule->n_points(); ++qp)
  {
    Real qp_temperature = 0;
    if (_has_temp)
    {
      qp_temperature = _temperature[qp];
      if (_temperature[qp] < 0)
      {
        std::stringstream msg;
        msg << "WARNING:  In HeatConductionMaterial3species:  negative temperature!\n"
            << "\tResetting to zero.\n"
            << "\t_qp: " << qp << "\n"
            << "\ttemp: " << _temperature[qp] << "\n"
            << "\telem: " << _current_elem->id() << "\n"
            << "\tproc: " << processor_id() << "\n";
        mooseWarning(msg.str());
        qp_temperature = 0;
      }
    }
    _thermal_conductivity[qp] = 0.0;
    _thermal_conductivity_dT[qp] = 0.0;
    if (_thermal_conductivity_temperature_function_1)
    {
      Point p;
      _thermal_conductivity[qp] += _mass_fraction_1[qp] *
          _thermal_conductivity_temperature_function_1->value(qp_temperature, p);
      _thermal_conductivity_dT[qp] += _mass_fraction_1[qp] *
          _thermal_conductivity_temperature_function_1->timeDerivative(qp_temperature, p);
    }
    else
    {
      _thermal_conductivity[qp] += _mass_fraction_1[qp] *
                                   _my_thermal_conductivity_1;
      _thermal_conductivity_dT[qp] += 0;
    }
    if (_thermal_conductivity_temperature_function_2)
    {
      Point p;
      _thermal_conductivity[qp] += _mass_fraction_2[qp] *
          _thermal_conductivity_temperature_function_2->value(qp_temperature, p);
      _thermal_conductivity_dT[qp] += _mass_fraction_2[qp] *
          _thermal_conductivity_temperature_function_2->timeDerivative(qp_temperature, p);
    }
    else
    {
      _thermal_conductivity[qp] += _mass_fraction_2[qp] *
                                   _my_thermal_conductivity_2;
      _thermal_conductivity_dT[qp] += 0;
    }
    if (_thermal_conductivity_temperature_function_3)
    {
      Point p;
      _thermal_conductivity[qp] += _mass_fraction_3[qp] *
          _thermal_conductivity_temperature_function_3->value(qp_temperature, p);
      _thermal_conductivity_dT[qp] += _mass_fraction_3[qp] *
          _thermal_conductivity_temperature_function_3->timeDerivative(qp_temperature, p);
    }
    else
    {
      _thermal_conductivity[qp] += _mass_fraction_3[qp] *
                                   _my_thermal_conductivity_3;
      _thermal_conductivity_dT[qp] += 0;
    }

    _specific_heat[qp] = 0.0;
    if (_specific_heat_temperature_function_1)
    {
      Point p;
      _specific_heat[qp] += _mass_fraction_1[qp] *
          _specific_heat_temperature_function_1->value(qp_temperature, p);
    }
    else
    {
      _specific_heat[qp] += _mass_fraction_1[qp] *
                            _my_specific_heat_1;
    }
    if (_specific_heat_temperature_function_2)
    {
      Point p;
      _specific_heat[qp] += _mass_fraction_2[qp] *
          _specific_heat_temperature_function_2->value(qp_temperature, p);
    }
    else
    {
      _specific_heat[qp] += _mass_fraction_2[qp] *
                            _my_specific_heat_2;
    }
    if (_specific_heat_temperature_function_3)
    {
      Point p;
      _specific_heat[qp] += _mass_fraction_3[qp] *
          _specific_heat_temperature_function_3->value(qp_temperature, p);
    }
    else
    {
      _specific_heat[qp] += _mass_fraction_3[qp] *
                            _my_specific_heat_3;
    }
  }
}
