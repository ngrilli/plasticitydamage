/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "LinearIsoElasticPFDamageMiehe.h"
#include "libmesh/utility.h"
#include "MathUtils.h"

template <>
InputParameters
validParams<LinearIsoElasticPFDamageMiehe>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Phase-field fracture model energy contribution to damage "
                             "growth-isotropic elasticity and undamaged stress under compressive "
                             "strain");
  params.addRequiredCoupledVar("b", "Laplacian of the order parameter for damage");
  params.addRequiredCoupledVar("initialc", "initial damage");
  params.addRequiredParam<Real>("l", "Interface width");
  params.addRequiredParam<Real>("visco", "Viscosity parameter");
  params.addRequiredParam<MaterialPropertyName>("gc_prop_var","Material property name with gc value");
  params.addParam<Real>("kdamage", 1e-6, "Stiffness of damaged matrix");

  return params;
}

LinearIsoElasticPFDamageMiehe::LinearIsoElasticPFDamageMiehe(const InputParameters & parameters)
  : ComputeStressBase(parameters),
    _gc_prop(getMaterialProperty<Real>("gc_prop_var")),
    _b(coupledValue("b")),
    _b_old(coupledValueOld("b")),
    _initialc(coupledValue("initialc")),
    _initialc_old(coupledValueOld("initialc")),
    _c(declareProperty<Real>("c")),
    _c_old(declarePropertyOld<Real>("c")),
    _damage(declareProperty<Real>("damage")),
    _damage_old(declarePropertyOld<Real>("damage")),
    _kdamage(getParam<Real>("kdamage")),
    _G0_pos(declareProperty<Real>("G0_pos")),
    _G0_pos_old(declarePropertyOld<Real>("G0_pos")),
    _x_bracket(declareProperty<Real>("x_bracket")),
    _x_bracket_old(declarePropertyOld<Real>("x_bracket")),
    _l(getParam<Real>("l")),
    _visco(getParam<Real>("visco")),
    _etens(LIBMESH_DIM),
    _epos(LIBMESH_DIM),
    _eigval(LIBMESH_DIM)
{
}

void
LinearIsoElasticPFDamageMiehe::computeQpStress()
{
  updateVar();
  updateJacobian();
}

void
LinearIsoElasticPFDamageMiehe::updateVar()
{
  // Isotropic elasticity is assumed
  Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  Real c_old = _initialc_old[_qp] + _damage_old[_qp];
  Real b_old = _b_old[_qp];
  const Real gc = _gc_prop[_qp];
  Real signx;
  Real xfac = _kdamage;
  if (c_old < 1.0)
    xfac += Utility::pow<2>(1.0 - c_old);

  _mechanical_strain[_qp].symmetricEigenvaluesEigenvectors(_eigval, _eigvec);

  // Tensors of outerproduct of eigen vectors
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    for (unsigned int j = 0; j < LIBMESH_DIM; ++j)
      for (unsigned int k = 0; k < LIBMESH_DIM; ++k)
        _etens[i](j, k) = _eigvec(j, i) * _eigvec(k, i);

  Real etr = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    etr += _eigval[i];

  Real etrpos = (std::abs(etr) + etr) / 2.0;
  Real etrneg = (std::abs(etr) - etr) / 2.0;

  RankTwoTensor stress0pos, stress0neg;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
  {
    stress0pos +=
        _etens[i] * (lambda * etrpos + 2.0 * mu * (std::abs(_eigval[i]) + _eigval[i]) / 2.0);
    stress0neg +=
        _etens[i] * (lambda * etrneg + 2.0 * mu * (std::abs(_eigval[i]) - _eigval[i]) / 2.0);
  }

  // Damage associated with positive component of stress
  _stress[_qp] = stress0pos * xfac - stress0neg;

  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    _epos[i] = (std::abs(_eigval[i]) + _eigval[i]) / 2.0;

  Real val = 0.0;
  for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    val += Utility::pow<2>(_epos[i]);
  val *= mu;

  // Energy with positive principal strains
  _G0_pos[_qp] = lambda * Utility::pow<2>(etrpos) / 2.0 + val;
  _x_bracket[_qp] = _l * _b_old[_qp] + 2.0 * (1.0 - c_old) * _G0_pos[_qp] / gc - c_old / _l;
  signx = MathUtils::sign(_x_bracket[_qp]);
  _x_bracket[_qp] = ((signx + 1.0) / 2.0) * _x_bracket[_qp] / _visco;
  _damage[_qp] = _damage_old[_qp] + _x_bracket[_qp] * _dt;

  if ( ( _damage[_qp] + _initialc_old[_qp] ) > 1.0 ) {
    _damage[_qp] = 1.0 - _initialc_old[_qp];
  }
  _c[_qp] = _damage[_qp] + _initialc_old[_qp];

}

void
LinearIsoElasticPFDamageMiehe::updateJacobian()
{
  _Jacobian_mult[_qp] = _elasticity_tensor[_qp];
}
