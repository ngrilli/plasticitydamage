/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef HEATCONDUCTIONMATERIAL3SPECIES_H
#define HEATCONDUCTIONMATERIAL3SPECIES_H

#include "Material.h"

// Forward Declarations
class HeatConductionMaterial3species;
class Function;

template <>
InputParameters validParams<HeatConductionMaterial3species>();

/**
 * Specific heat and thermal conductivity for a material
 * with three different chemical species
 */
class HeatConductionMaterial3species : public Material
{
public:
  HeatConductionMaterial3species(const InputParameters & parameters);

protected:
  virtual void computeProperties();

  const bool _has_temp;
  const VariableValue & _temperature;

// mass fractions of reactants
  const bool _has_mass_fraction_1;
  const VariableValue & _mass_fraction_1;
  const bool _has_mass_fraction_2;
  const VariableValue & _mass_fraction_2;
  const bool _has_mass_fraction_3;
  const VariableValue & _mass_fraction_3;

  const Real _my_thermal_conductivity_1;
  const Real _my_specific_heat_1;
  const Real _my_thermal_conductivity_2;
  const Real _my_specific_heat_2;
  const Real _my_thermal_conductivity_3;
  const Real _my_specific_heat_3;

  MaterialProperty<Real> & _thermal_conductivity;
  MaterialProperty<Real> & _thermal_conductivity_dT;
  Function * _thermal_conductivity_temperature_function_1;
  Function * _thermal_conductivity_temperature_function_2;
  Function * _thermal_conductivity_temperature_function_3;

  MaterialProperty<Real> & _specific_heat;
  Function * _specific_heat_temperature_function_1;
  Function * _specific_heat_temperature_function_2;
  Function * _specific_heat_temperature_function_3;
};

#endif // HEATCONDUCTIONMATERIAL3SPECIES_H
