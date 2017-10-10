/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRHENIUSHEATENERGYRATELIMIT_H
#define ARRHENIUSHEATENERGYRATELIMIT_H

#include "DerivativeMaterialInterface.h"
#include "Kernel.h"

// Forward Declarations
class ArrheniusHeatEnergyRateLimit;

template <>
InputParameters validParams<ArrheniusHeatEnergyRateLimit>();

/**
 * Provides a heat source from chemical reactions:
 * Arrhenius exponential law
 * upper threshold for the reaction rate
 */
class ArrheniusHeatEnergyRateLimit : public DerivativeMaterialInterface<Kernel>
{
public:
  ArrheniusHeatEnergyRateLimit(const InputParameters & parameters);

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

  /// d(mass_fraction_rate)/d(T)
  const MaterialProperty<Real> & _dmass_fraction_rate_dtemperature;

  /// d(mass_fraction_rate)/d(mass_fraction)
  const MaterialProperty<Real> & _dmass_fraction_rate_dmass_fraction;

  /// Heat of combustion (Q)
  Real _combustion_heat;

  /// Density
  const MaterialProperty<Real> & _density;
};

#endif // ARRHENIUSHEATENERGYRATELIMIT_H
