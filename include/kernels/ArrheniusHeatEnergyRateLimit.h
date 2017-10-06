/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRHENIUSHEATENERGYRATELIMIT_H
#define ARRHENIUSHEATENERGYRATELIMIT_H

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
class ArrheniusHeatEnergyRateLimit : public Kernel
{
public:
  ArrheniusHeatEnergyRateLimit(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// optional parameter that allows multiple mechanics models to be defined
  std::string _base_name;

  /// Auxiliary variable: mass fraction of reactant
  const VariableValue & _mass_fraction;
  const unsigned int _mass_fraction_var;

  Real _exponential_prefactor;
  Real _exponential_coefficient;
  Real _exponential_factor;

  /// Upper threshold for the reaction rate
  Real _rate_limit;
};

#endif // ARRHENIUSHEATENERGYRATELIMIT_H
