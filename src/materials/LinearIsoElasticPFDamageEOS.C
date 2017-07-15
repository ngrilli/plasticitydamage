/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LinearIsoElasticPFDamageEOS.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<LinearIsoElasticPFDamageEOS>()
{
  InputParameters params = validParams<LinearIsoElasticPFDamage>();
  params.addClassDescription("Phase-field fracture model energy contribution to damage growth-isotropic elasticity and undamaged stress under compressive strain and Birch-Murnaghan EOS");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");

  return params;
}

LinearIsoElasticPFDamageEOS::LinearIsoElasticPFDamageEOS(const InputParameters & parameters) :
    LinearIsoElasticPFDamage(parameters),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref"))
{
}

void
LinearIsoElasticPFDamageEOS::updateVar()
{
  RankTwoTensor stress0pos, stress0neg, stress0;
  //Isotropic elasticity is assumed
  Real lambda = _elasticity_tensor[_qp](0,0,1,1);
  Real mu = _elasticity_tensor[_qp](0,1,0,1);
  Real vv0 = 1.0 + _mechanical_strain[_qp].trace(); //specific volume
  Real c = _c[_qp];
  Real xfac = Utility::pow<2>(1.0-c) + _kdamage;

  //EOS
  lambda = lambda / std::pow( vv0 , _n_Murnaghan + 1.0 );
  mu = mu / std::pow( vv0 , _n_Murnaghan + 1.0 );

  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

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
  //_stress[_qp].addIa(-1.0/3.0 * _stress[_qp].trace());
  //_stress[_qp].addIa( (_Bulk_Modulus_Ref/_n_Murnaghan) * ( 1.0 - std::pow( 1.0/(1.0+etr) , _n_Murnaghan ) ) );

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
