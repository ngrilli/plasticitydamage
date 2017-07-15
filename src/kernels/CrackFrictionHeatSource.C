#include "CrackFrictionHeatSource.h"

template <>
InputParameters
validParams<CrackFrictionHeatSource>()
{
  InputParameters params = validParams<HeatSource>();
  params.addClassDescription("Crack friction heat source kernel");
  //params.addCoupledVar("elec", "Electric potential for joule heating.");
  //params.addCoupledVar("args", "Vector of arguments of the diffusivity");
  params.addRequiredParam<Real>("friction_coefficient","crack friction coefficient");
  params.addRequiredCoupledVar("dcdx","First derivative of damage with respect to x");
  params.addRequiredCoupledVar("dcdy","First derivative of damage with respect to y");
  //params.addRequiredCoupledVar("dcdz","First derivative of damage with respect to z");
  return params;
}

CrackFrictionHeatSource::CrackFrictionHeatSource(const InputParameters & parameters)
  : HeatSource(parameters),
    //_grad_elec(coupledGradient("elec")),
    //_elec_var(coupled("elec")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _dcdx(coupledValue("dcdx")),
    _dcdy(coupledValue("dcdy")),
    //_dcdz(coupledValue("dcdz")),
    _crack_normal(getMaterialProperty<std::vector<Real>>("crack_normal")), // Normal to the crack surface
    _crack_normal_norm(getMaterialProperty<Real>("crack_normal_norm")),
    _stress(getMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _strain_rate(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_rate")),
    _Jacobian_mult(getMaterialProperty<RankFourTensor>(_base_name + "Jacobian_mult")),
    _friction_coefficient(getParam<Real>("friction_coefficient")),
    _friction_force(getMaterialProperty<std::vector<Real>>("friction_force")),
    _slide_velocity(getMaterialProperty<std::vector<Real>>("slide_velocity")),
    _heat_source_rate(getMaterialProperty<Real>("heat_source_rate")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp)
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);
}

Real
CrackFrictionHeatSource::computeQpResidual()
{
  return -_heat_source_rate[_qp] * _test[_i][_qp];
}

Real
CrackFrictionHeatSource::computeQpJacobian()
{
  return -0.0 * _phi[_j][_qp] * _test[_i][_qp];
}

Real
CrackFrictionHeatSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  Real val;
  unsigned int i, j, l, h;

  val = 0.0;
  if ( _heat_source_rate[_qp] > 0.0 ) {
    for (unsigned int k = 0; k < _ndisp; ++k)
    {
      if (jvar == _disp_var[k]) {
        for (unsigned int h = 0; h < 2; ++h) {
          for (unsigned int j = 0; j < 2; ++j) {
            for (unsigned int i = 0; i < 2; ++i) {
              val += _Jacobian_mult[_qp](i,j,k,h) * _grad_phi[_j][_qp](h) * _crack_normal[_qp][j] * _slide_velocity[_qp][i];
            }
          }
        }
        for (unsigned int h = 0; h < 2; ++h) {
          val += _friction_force[_qp][k] * _grad_phi[_j][_qp](h) * _crack_normal[_qp][h] / ( 2.0 * _dt );
        }
        for (unsigned int h = 0; h < 2; ++h) {
          val += _friction_force[_qp][h] * _grad_phi[_j][_qp](h) * _crack_normal[_qp][k] / ( 2.0 * _dt );
        }
        val = -_friction_coefficient * val * _test[_i][_qp];
      }
    }
  }
  return val;
}
