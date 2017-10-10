/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRHENIUSMASSFRACTIONRATELIMIT_H
#define ARRHENIUSMASSFRACTIONRATELIMIT_H

#include "DerivativeMaterialInterface.h"
#include "Kernel.h"

// Forward Declarations
class ArrheniusMassFractionRateLimit;

template <>
InputParameters validParams<ArrheniusMassFractionRateLimit>();

/**
 * Provides the rate of reactant
 * mass fraction decrease during chemical reactions:
 * Arrhenius exponential law
 * with an upper threshold of the reaction rate
 */
class ArrheniusMassFractionRateLimit : public DerivativeMaterialInterface<Kernel>
{
public:
  ArrheniusMassFractionRateLimit(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// optional parameter that allows multiple mechanics models to be defined
  std::string _base_name;

  /// Auxiliary variable: mass_fraction
  const VariableValue & _mass_fraction;

  /// MOOSE variable number for the mass fraction variable
  const unsigned int _mass_fraction_var;

  /// name of mass fraction variable
  VariableName _mass_fraction_name;

  /// Auxiliary variable: temperature
  const VariableValue & _temperature;

  /// MOOSE variable number for the temperature variable
  const unsigned int _temperature_var;

  /// name of temperature variable
  VariableName _temperature_name;

  /// Decrease rate of the mass fraction due to chemistry
  const MaterialProperty<Real> & _mass_fraction_rate;

  /// Derivative of the decrease rate of the mass fraction due to chemistry with respect to temperature
  const MaterialProperty<Real> & _dmass_fraction_rate_dtemperature;

  /// Derivative of the decrease rate of the mass fraction due to chemistry with respect to mass fraction
  const MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction;
};

#endif // ARRHENIUSMASSFRACTIONRATELIMIT_H
