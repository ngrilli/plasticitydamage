/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRHENIUSMASSFRACTIONRATELIMIT_H
#define ARRHENIUSMASSFRACTIONRATELIMIT_H

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
class ArrheniusMassFractionRateLimit : public Kernel
{
public:
  ArrheniusMassFractionRateLimit(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// optional parameter that allows multiple mechanics models to be defined
  std::string _base_name;

  /// Auxiliary variable: temperature
  const VariableValue & _temperature;
  const unsigned int _temperature_var;

  Real _exponential_prefactor;
  Real _exponential_coefficient;
  Real _exponential_factor;

  /// Upper threshold for the reaction rate
  Real _rate_limit;
};

#endif // ARRHENIUSMASSFRACTIONRATELIMIT_H
