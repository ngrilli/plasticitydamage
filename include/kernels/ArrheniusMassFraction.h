/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRHENIUSMASSFRACTION_H
#define ARRHENIUSMASSFRACTION_H

#include "Kernel.h"

// Forward Declarations
class ArrheniusMassFraction;

template <>
InputParameters validParams<ArrheniusMassFraction>();

/**
 * Provides the rate of reactant
 * mass fraction decrease during chemical reactions:
 * Arrhenius exponential law
 */
class ArrheniusMassFraction : public Kernel
{
public:
  ArrheniusMassFraction(const InputParameters & parameters);

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
};

#endif // ARRHENIUSMASSFRACTION_H
