/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef ARRHENIUSHEATENERGY_H
#define ARRHENIUSHEATENERGY_H

#include "Kernel.h"

// Forward Declarations
class ArrheniusHeatEnergy;

template <>
InputParameters validParams<ArrheniusHeatEnergy>();

/**
 * Provides a heat source from chemical reactions:
 * Arrhenius exponential law
 */
class ArrheniusHeatEnergy : public Kernel
{
public:
  ArrheniusHeatEnergy(const InputParameters & parameters);

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
};

#endif // ARRHENIUSHEATENERGY_H
