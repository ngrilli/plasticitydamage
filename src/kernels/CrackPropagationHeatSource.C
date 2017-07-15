#include "CrackPropagationHeatSource.h"

template <>
InputParameters
validParams<CrackPropagationHeatSource>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Crack propagation heat source kernel");
  params.addRequiredCoupledVar("c","Damage");
  params.addRequiredParam<MaterialPropertyName>(
      "G0_var", "Material property name with undamaged strain energy driving damage (G0_pos)");
  return params;
}

CrackPropagationHeatSource::CrackPropagationHeatSource(const InputParameters & parameters)
  : HeatSource(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _c(coupledValue("c")),
    _c_old(coupledValueOld("c")),
    _c_var(coupled("c")),
    _G0_pos(getMaterialProperty<Real>("G0_var")),
    _dG0_pos_dstrain(isParamValid("dG0_dstrain_var")
                         ? &getMaterialProperty<RankTwoTensor>("dG0_dstrain_var")
                         : NULL),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
CrackPropagationHeatSource::computeQpResidual()
{
  Real ddot, heat_source_rate;

  ddot = ( _c[_qp] - _c_old[_qp] ) / _dt;
  heat_source_rate = 2.0 * ( 1.0 - _c[_qp] ) * _G0_pos[_qp] * ddot;

  if ( heat_source_rate < 0.0 ) {
    return 0.0;
  }
  else {
    return - heat_source_rate * _test[_i][_qp];
  }
}

Real
CrackPropagationHeatSource::computeQpJacobian()
{
  return -0.0 * _phi[_j][_qp] * _test[_i][_qp];
}

Real
CrackPropagationHeatSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val, ddot, dRdd;
  unsigned int i, j, l, h, c_comp;
  bool disp_flag = false;

  ddot = ( _c[_qp] - _c_old[_qp] ) / _dt;
  dRdd = 2.0 * _G0_pos[_qp] * ( ddot - ( 1.0 - _c[_qp] ) / _dt );

  val = 0.0;
  // Contribution of auxiliary variable to off diag Jacobian of temperature
  for (unsigned int k = 0; k < _ndisp; ++k)
  {
    if (jvar == _c_var)
      return dRdd * _phi[_j][_qp] * _test[_i][_qp];
    else if (jvar == _disp_var[k])
    {
      c_comp = k;
      disp_flag = true;
    }
  }

  // Contribution of displacements to off diag Jacobian of c
  if (disp_flag && _dG0_pos_dstrain != NULL)
  {
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      val += ((*_dG0_pos_dstrain)[_qp](c_comp, i) + (*_dG0_pos_dstrain)[_qp](i, c_comp)) * _grad_phi[_j][_qp](i);
    }
    val *= - ( 1.0 - _c[_qp] ) * ddot;
  }

  return val;

}
