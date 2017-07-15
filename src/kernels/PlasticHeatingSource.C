#include "PlasticHeatingSource.h"

template <>
InputParameters
validParams<PlasticHeatingSource>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Heat source kernel from plasticity");
  params.addParam<MaterialPropertyName>("W0p","plastic energy");
  return params;
}

PlasticHeatingSource::PlasticHeatingSource(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _W0p(getMaterialProperty<Real>("W0p")),
    _W0p_old(getMaterialPropertyOld<Real>("W0p")),
    _dW0p_dstrain(isParamValid("dW0p_dstrain_var") ? &getMaterialProperty<RankTwoTensor>("dW0p_dstrain_var"): NULL),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
PlasticHeatingSource::computeQpResidual()
{
  return - ( _W0p[_qp] - _W0p_old[_qp] ) / _dt * _test[_i][_qp];
}

Real
PlasticHeatingSource::computeQpJacobian()
{
  return 0.0 * _phi[_j][_qp] * _test[_i][_qp];
}

Real
PlasticHeatingSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  unsigned int c_comp;
  bool disp_flag = false;

  val = 0.0;
  // Contribution of auxiliary variable to off diag Jacobian of temperature
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
      val += ( (*_dW0p_dstrain)[_qp](c_comp, i) + (*_dW0p_dstrain)[_qp](i, c_comp) ) * _grad_phi[_j][_qp](i) / _dt;
    }
  }

  return val;
}
