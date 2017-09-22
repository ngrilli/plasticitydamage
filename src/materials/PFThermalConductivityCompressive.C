/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PFThermalConductivityCompressive.h"

// libmesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters validParams<PFThermalConductivityCompressive>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("This material calculates effective thermal "
                             "conductivity of partially damaged material "
                             "k_eff = (1-c^2)*k_m + c^2*k_c, "
                             "the conductivity is stress dependent, "
                             "if the material is under compression, "
                             "then the conductivity is not degraded");
  params.addRequiredCoupledVar("c", "Coupled phase field");
  params.addRequiredParam<Real>("material_thermal_conductivity",
                                "Thermal conductivity of the undamaged material");
  params.addParam<Real>("crack_thermal_conductivity", 0.0,
                        "Thermal conductivity of fully damaged material");
  return params;
}

PFThermalConductivityCompressive::PFThermalConductivityCompressive(const InputParameters & parameters)
  : Material(parameters),
    _c(coupledValue("c")),
    _k_m(getParam<Real>("material_thermal_conductivity")),
    _k_c(getParam<Real>("crack_thermal_conductivity")),
    _k(declareProperty<Real>("thermal_conductivity")),
    _stress(getMaterialPropertyByName<RankTwoTensor>("stress")),
    _friction_normal_force(getMaterialPropertyByName<Real>("friction_normal_force"))
{
}

void PFThermalConductivityCompressive::computeQpProperties()
{
    Real qp_c = _c[_qp];

    if (_c[_qp] < 0.0) qp_c = 0.0;
    if (_c[_qp] > 1.0) qp_c = 1.0;

    if ( _friction_normal_force[_qp] <= 0.0 ) { // compressive load
      _k[_qp] = _k_m;
    }
    else { // tensile load
      _k[_qp] = _k_m * (1.0 - qp_c * qp_c) + _k_c * qp_c * qp_c;
    }
}
