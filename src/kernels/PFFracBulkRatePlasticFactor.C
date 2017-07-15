/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "PFFracBulkRatePlasticFactor.h"
#include "MathUtils.h"

template<>
InputParameters validParams<PFFracBulkRatePlasticFactor>()
{
  InputParameters params = validParams<KernelValue>();
  params.addClassDescription("Kernel to compute bulk energy contribution to damage order parameter residual equation with plasticity");
  params.addRequiredParam<Real>("l","Interface width");
  params.addRequiredParam<Real>("visco","Viscosity parameter");
  params.addRequiredParam<Real>("Wc","Threshold fracture energy");
  params.addRequiredParam<Real>("plastic_factor","Prefactor of the plastic contribution to damage");
  params.addRequiredParam<MaterialPropertyName>("gc_prop_var", "Material property name with gc value");
  params.addRequiredParam<MaterialPropertyName>("W0e_var", "Elastic work");
  params.addRequiredParam<MaterialPropertyName>("W0p_var", "Plastic work");
  params.addParam<MaterialPropertyName>("dW0e_dstrain_var", "Derivative of W0e with strain");
  params.addParam<MaterialPropertyName>("dW0p_dstrain_var", "Derivative of W0e with strain");
  params.addRequiredCoupledVar("beta", "Auxiliary variable");
  params.addCoupledVar("displacements",
                       "The string of displacements suitable for the problem statement");
  params.addParam<std::string>("base_name", "Material property base name");

  return params;
}

PFFracBulkRatePlasticFactor::PFFracBulkRatePlasticFactor(const InputParameters & parameters):
  KernelValue(parameters),
  _gc_prop(getMaterialProperty<Real>("gc_prop_var")),
  _W0e(getMaterialProperty<Real>("W0e_var")),
  _W0p(getMaterialProperty<Real>("W0p_var")),
  _dW0e_dstrain(isParamValid("dW0e_dstrain_var") ? &getMaterialProperty<RankTwoTensor>("dW0e_dstrain_var"): NULL),
  _dW0p_dstrain(isParamValid("dW0p_dstrain_var") ? &getMaterialProperty<RankTwoTensor>("dW0p_dstrain_var"): NULL),
  _betaval(coupledValue("beta")),
  _beta_var(coupled("beta")),
  _ndisp(coupledComponents("displacements")),
  _disp_var(_ndisp),
  _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
  _l(getParam<Real>("l")),
  _visco(getParam<Real>("visco")),
  _Wc(getParam<Real>("Wc")),
  _plastic_factor(getParam<Real>("plastic_factor"))
{
}

Real
PFFracBulkRatePlasticFactor::precomputeQpResidual()
{
  const Real gc = _gc_prop[_qp];
  const Real c = _u[_qp];
  const Real x = _l * _betaval[_qp] + 2.0 * (1.0 - c) * (_W0e[_qp]+_plastic_factor * _W0p[_qp]-_Wc) / gc - c / _l;

  return -((std::abs(x) + x) / 2.0) / _visco;
}

Real
PFFracBulkRatePlasticFactor::precomputeQpJacobian()
{
  const Real gc = _gc_prop[_qp];
  const Real c = _u[_qp];
  const Real x = _l * _betaval[_qp] + 2.0 * (1.0 - c) * (_W0e[_qp]+_plastic_factor * _W0p[_qp]-_Wc) / gc - c / _l;

  return (MathUtils::sign(x) + 1.0) / 2.0 * (2.0 * (_W0e[_qp]+_plastic_factor * _W0p[_qp]-_Wc) / gc + 1.0 / _l) / _visco *
         _phi[_j][_qp];
}

Real
PFFracBulkRatePlasticFactor::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int c_comp;
  bool disp_flag = false;

  const Real c = _u[_qp];
  const Real gc = _gc_prop[_qp];

  const Real x = _l * _betaval[_qp] + 2.0 * (1.0 - c) * ( (_W0e[_qp]+_plastic_factor * _W0p[_qp]-_Wc) / gc) - c / _l;

  const Real signx = MathUtils::sign(x);

  Real xfacbeta = -((signx + 1.0) / 2.0) / _visco * _l;
  Real xfac = -((signx + 1.0) / 2.0) / _visco * 2.0 * (1.0 - c) / gc;

  // Contribution of auxiliary variable to off diag Jacobian of c
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _beta_var)
      return xfacbeta * _phi[_j][_qp] * _test[_i][_qp];
    else if (jvar == _disp_var[k])
    {
      c_comp = k;
      disp_flag = true;
    }
  }

  // Contribution of displacements to off diag Jacobian of c
  if (disp_flag && _dW0e_dstrain != NULL && _dW0p_dstrain != NULL)
  {
    Real val = 0.0;
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      val += ((*_dW0e_dstrain)[_qp](c_comp, i) + (*_dW0e_dstrain)[_qp](i, c_comp)) / 2.0 *
             _grad_phi[_j][_qp](i);
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
      val += _plastic_factor * ((*_dW0p_dstrain)[_qp](c_comp, i) + (*_dW0p_dstrain)[_qp](i, c_comp)) / 2.0 *
             _grad_phi[_j][_qp](i);

    return xfac * val * _test[_i][_qp];
  }

  return 0.0;
}
