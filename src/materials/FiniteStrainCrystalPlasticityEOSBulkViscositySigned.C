/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainCrystalPlasticityEOSBulkViscositySigned.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityEOSBulkViscositySigned>()
{
  InputParameters params = validParams<FiniteStrainCrystalPlasticity>();
  params.addClassDescription("Crystal Plasticity base class: FCC system with power law flow rule implemented: bulk viscosity and equation of state");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");

  return params;
}

FiniteStrainCrystalPlasticityEOSBulkViscositySigned::FiniteStrainCrystalPlasticityEOSBulkViscositySigned(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _temp(coupledValue("temp")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature"))
{
}

void
FiniteStrainCrystalPlasticityEOSBulkViscositySigned::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  Real trD;

  iden.zero();
  iden.addIa(1.0);

  _fe = _dfgrd_tmp * _fp_prev_inv; // _fp_inv  ==> _fp_prev_inv

  ce = _fe.transpose() * _fe;
  ce_pk2 = ce * _pk2_tmp;
  ce_pk2 = ce_pk2 / _fe.det();

  // Calculate Schmid tensor and resolved shear stresses
  for (unsigned int i = 0; i < _nss; ++i)
    _tau(i) = ce_pk2.doubleContraction(_s0[i]);

  getSlipIncrements(); // Calculate dslip,dslipdtau

  if (_err_tol)
    return;

  eqv_slip_incr.zero();
  for (unsigned int i = 0; i < _nss; ++i)
    eqv_slip_incr += _s0[i] * _slip_incr(i);

  eqv_slip_incr = iden - eqv_slip_incr;
  _fp_inv = _fp_old_inv * eqv_slip_incr;
  _fe = _dfgrd_tmp * _fp_inv;

  ce = _fe.transpose() * _fe;
  ee = ce - iden;
  ee *= 0.5;

  pk2_new = _elasticity_tensor[_qp] * ee;

  pk2_new.addIa(-1.0/3.0 * pk2_new.trace());
  pk2_new.addIa( - (_Bulk_Modulus_Ref/_n_Murnaghan) * std::pow( 1.0/_fe.det() , _n_Murnaghan ) );
  pk2_new.addIa( (_Bulk_Modulus_Ref/_n_Murnaghan) * std::exp( - _n_Murnaghan * _thermal_expansion * ( _temp[_qp] - _reference_temperature ) ) );

  trD = ( _deformation_gradient[_qp].det() - _deformation_gradient_old[_qp].det() ) / _dt;
  trD /= _deformation_gradient_old[_qp].det();

  pk2_new.addIa( _C0 * trD * abs(trD) );
  pk2_new.addIa( _C1 * trD );

  resid = _pk2_tmp - pk2_new;
}
