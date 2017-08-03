#include "PlasticHeatingSourceMiehe2016.h"

template <>
InputParameters
validParams<PlasticHeatingSourceMiehe2016>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Heat source kernel from plasticity in Miehe 2016");
  params.addRequiredCoupledVar("c", "phase field damage");
  params.addParam<MaterialPropertyName>("W0p","plastic energy");
  return params;
}

PlasticHeatingSourceMiehe2016::PlasticHeatingSourceMiehe2016(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _W0p(getMaterialProperty<Real>("W0p")),
    _W0p_old(getMaterialPropertyOld<Real>("W0p")),
    _dW0p_dstrain(isParamValid("dW0p_dstrain_var") ? &getMaterialProperty<RankTwoTensor>("dW0p_dstrain_var"): NULL),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _cval(coupledValue("c")),
    _c_var(coupled("c"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
PlasticHeatingSourceMiehe2016::computeQpResidual()
{
  return - std::pow( 1.0 - _cval[_qp] , 2.0) * ( _W0p[_qp] - _W0p_old[_qp] ) / _dt * _test[_i][_qp];
}

Real
PlasticHeatingSourceMiehe2016::computeQpJacobian()
{
  return 0.0 * _phi[_j][_qp] * _test[_i][_qp];
}

Real
PlasticHeatingSourceMiehe2016::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  unsigned int c_comp;
  bool disp_flag = false;

  val = 0.0;
  // Contribution of auxiliary variable to off diag Jacobian of temperature
  if (jvar == _c_var)
  {
    return 2.0 * (1.0 - _cval[_qp]) * ( _W0p[_qp] - _W0p_old[_qp] ) / _dt * _phi[_j][_qp] * _test[_i][_qp];
  }
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k])
    {
      c_comp = k;
      disp_flag = true;
    }
  }

  // Contribution of displacements to off diag Jacobian of c
  if (disp_flag && _dW0p_dstrain != NULL)
  {
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      val += ( (*_dW0p_dstrain)[_qp](c_comp, i) + (*_dW0p_dstrain)[_qp](i, c_comp) ) * _grad_phi[_j][_qp](i) / _dt / 2.0;
    }
  }

  return - val * std::pow( 1.0 - _cval[_qp] , 2.0) * _test[_i][_qp];
}
