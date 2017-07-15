/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LinearIsoElasticPFDamageEOSBulkViscosity.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<LinearIsoElasticPFDamageEOSBulkViscosity>()
{
  InputParameters params = validParams<LinearIsoElasticPFDamage>();
  params.addClassDescription("Phase-field fracture model energy contribution to damage growth-isotropic elasticity and undamaged stress under compressive strain and Birch-Murnaghan EOS and bulk viscosity");
  params.addRequiredCoupledVar("temp","Temperature");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("C0", "Von Neumann coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");
  params.addRequiredParam<Real>("thermal_expansion", "Thermal expansion coefficient");
  params.addRequiredParam<Real>("reference_temperature", "reference temperature for thermal expansion");

  return params;
}

LinearIsoElasticPFDamageEOSBulkViscosity::LinearIsoElasticPFDamageEOSBulkViscosity(const InputParameters & parameters) :
    LinearIsoElasticPFDamage(parameters),
    _temp(coupledValue("temp")),
    _strain_rate(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_rate")),
    _mechanical_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _mechanical_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "mechanical_strain")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _thermal_expansion(getParam<Real>("thermal_expansion")),
    _reference_temperature(getParam<Real>("reference_temperature"))
{
}

void
LinearIsoElasticPFDamageEOSBulkViscosity::updateVar()
{
  RankTwoTensor stress0pos, stress0neg;
  RankTwoTensor FSmallStrain;
  //Isotropic elasticity is assumed
  Real lambda = _elasticity_tensor[_qp](0,0,1,1);
  Real mu = _elasticity_tensor[_qp](0,1,0,1);
  Real c = _c[_qp];
  Real xfac = _kdamage;
  if (c < 1.0)
    xfac += Utility::pow<2>(1.0 - c);
  Real trD, specificVolume, Jdot;

  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

  //_strain_rate[_qp] = ( _mechanical_strain[_qp] - _mechanical_strain_old[_qp] ) / _dt;

  specificVolume = 1.0 + _mechanical_strain[_qp].trace();

  Jdot = _strain_rate[_qp].trace();

  trD = Jdot / specificVolume;

  FSmallStrain.zero();
  FSmallStrain.addIa( 1.0 );
  FSmallStrain += _mechanical_strain[_qp];

  //Tensors of outerproduct of eigen vectors
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _etens[i](j,k) = _eigvec(j,i) * _eigvec(k,i);

  Real etr = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    etr += _eigval[i];

  Real etrpos = (std::abs(etr) + etr) / 2.0;
  Real etrneg = (std::abs(etr) - etr) / 2.0;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    stress0pos += _etens[i] * (lambda * etrpos + 2.0 * mu * (std::abs(_eigval[i]) + _eigval[i]) / 2.0);
    stress0neg += _etens[i] * (lambda * etrneg + 2.0 * mu * (std::abs(_eigval[i]) - _eigval[i]) / 2.0);
  }

  //Damage associated with positive component of stress
  _stress[_qp] = stress0pos * xfac - stress0neg;

  //EOS
  _stress[_qp].addIa(-1.0/3.0 * _stress[_qp].trace());
  _stress[_qp].addIa( - (_Bulk_Modulus_Ref/_n_Murnaghan) * std::pow( 1.0/specificVolume , _n_Murnaghan ) );
  _stress[_qp].addIa( (_Bulk_Modulus_Ref/_n_Murnaghan) * std::exp( - _n_Murnaghan * _thermal_expansion * ( _temp[_qp] - _reference_temperature ) ) );

  //Bulk viscosity
  _stress[_qp].addIa( _C0 * trD * abs(trD) );
  _stress[_qp].addIa( _C1 * trD );

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    _epos[i] = (std::abs(_eigval[i]) + _eigval[i]) / 2.0;

  Real val = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    val += Utility::pow<2>(_epos[i]);
  val *= mu;

  //Energy with positive principal strains
  _G0_pos[_qp] = lambda * Utility::pow<2>(etrpos) / 2.0 + val;
  //Used in PFFracBulkRate Jacobian
  _dG0_pos_dstrain[_qp] = stress0pos;
  //Used in StressDivergencePFFracTensors Jacobian
  _dstress_dc[_qp] = -stress0pos * (2.0 * (1.0 - c));
}
