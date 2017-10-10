/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEARRHENIUSMASSFRACTIONRATELIMIT_H
#define COMPUTEARRHENIUSMASSFRACTIONRATELIMIT_H

#include "DerivativeMaterialInterface.h"
#include "Material.h"

/**
 * ComputeArrheniusMassFractionRateLimit computes the mass fraction rate
 * given by chemical reactions, upper threshold rate is used
 */
class ComputeArrheniusMassFractionRateLimit : public DerivativeMaterialInterface<Material>
{
public:
  ComputeArrheniusMassFractionRateLimit(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// optional parameter that allows multiple mechanics materials to be defined
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

  Real _exponential_prefactor;
  Real _exponential_coefficient;
  Real _exponential_factor;

  /// Upper threshold for the reaction rate
  Real _rate_limit;

  /// computed property: mass fraction rate
  MaterialProperty<Real> & _mass_fraction_rate;

  /// d(mass_fraction_rate)/d(T)
  MaterialProperty<Real> & _dmass_fraction_rate_dtemperature;

  /// d(mass_fraction_rate)/d(mass_fraction)
  MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction;
};

#endif // COMPUTEARRHENIUSMASSFRACTIONRATELIMIT_H
