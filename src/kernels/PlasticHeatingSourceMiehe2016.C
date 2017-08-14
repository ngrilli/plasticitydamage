#include "PlasticHeatingSourceMiehe2016.h"

template <>
InputParameters
validParams<PlasticHeatingSourceMiehe2016>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Heat source kernel from plasticity in Miehe 2016");
  params.addParam<MaterialPropertyName>("stress","Cauchy stress");
  params.addParam<MaterialPropertyName>("deformation_gradient","deformation gradient");
  params.addParam<MaterialPropertyName>("fp","plastic deformation gradient");
  return params;
}

PlasticHeatingSourceMiehe2016::PlasticHeatingSourceMiehe2016(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getMaterialProperty<RankTwoTensor>("stress")),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>("deformation_gradient")),
    _fp(getMaterialProperty<RankTwoTensor>("fp")),
    _fp_old(getMaterialPropertyOld<RankTwoTensor>("fp")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
PlasticHeatingSourceMiehe2016::computeQpResidual()
{
  Real detFe;
  RankTwoTensor WpToTrace, _fp_inv, _fe, _fe_inv;

  // invert plastic deformation gradient
  _fp_inv = _fp[_qp].inverse();

  // calculate elastic deformation gradient
  _fe = _deformation_gradient[_qp] * _fp_inv;
  _fe_inv = _fe.inverse();
  detFe = _fe.det();

  // calculate plastic work (energy per unit volume)
  WpToTrace = _stress[_qp] * _fe * ( _fp[_qp] - _fp_old[_qp] ) * _fp_inv * _fe_inv * detFe;

  // plastic work rate
  return -WpToTrace.trace() / _dt * _test[_i][_qp];
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
  //if (jvar == _c_var) _stress is actually coupled with _c, but for the moment we neglect
  //{
  //  return 2.0 * (1.0 - _cval[_qp]) * ( _W0p[_qp] - _W0p_old[_qp] ) / _dt * _phi[_j][_qp] * _test[_i][_qp];
  //}
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _disp_var[k])
    {
      c_comp = k;
      disp_flag = true;
    }
  }

  // Contribution of displacements to off diag Jacobian of the plastic work rate
  if (disp_flag)
  {
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      val += _stress[_qp](c_comp, i) * _grad_phi[_j][_qp](i) / _dt;
    }
  }

  return -val * _test[_i][_qp];
}
