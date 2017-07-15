/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "FiniteStrainCrystalPlasticityMcAuliffe.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"

template<>
InputParameters validParams<FiniteStrainCrystalPlasticityMcAuliffe>()
{
  InputParameters params = validParams<FiniteStrainCrystalPlasticity>();
  params.addClassDescription("Crystal Plasticity class: McAuliffe and Waisman 2015");
  params.addRequiredCoupledVar("c","Order parameter for damage");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");

  return params;
}

FiniteStrainCrystalPlasticityMcAuliffe::FiniteStrainCrystalPlasticityMcAuliffe(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
    _c(coupledValue("c")),
    _kdamage(getParam<Real>("kdamage")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _W0e(declareProperty<Real>("W0e")), // elastic energy
    _W0p(declareProperty<Real>("W0p")), // plastic energy
    _W0p_old(declarePropertyOld<Real>("W0p")), // plastic energy of previous increment
    _dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dW0e_dstrain(declareProperty<RankTwoTensor>("dW0e_dstrain")),
    _dW0p_dstrain(declareProperty<RankTwoTensor>("dW0p_dstrain")),
    _pk2_undamaged(declareProperty<RankTwoTensor>("pk2_undamaged")) // undamaged 2nd Piola Kirchoff Stress
{
}

void
FiniteStrainCrystalPlasticityMcAuliffe::preSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss_tmp = _gss_old[_qp];
    _W0p_tmp = _W0p_old[_qp];
    _accslip_tmp_old = _acc_slip_old[_qp];
  }
  else
  {
    if (_first_step_iter)
    {
      _gss_tmp = _gss_tmp_old = _gss_old[_qp];
      _W0p_tmp = _W0p_tmp_old = _W0p_old[_qp];
      _accslip_tmp_old = _acc_slip_old[_qp];
    }
    else
    {
      _gss_tmp = _gss_tmp_old;
      _W0p_tmp = _W0p_tmp_old;
    }
  }
}

void
FiniteStrainCrystalPlasticityMcAuliffe::solveStatevar()
{
  Real gmax, gdiff;
  unsigned int iterg;
  std::vector<Real> gss_prev(_nss);

  gmax = 1.1 * _gtol;
  iterg = 0;

  while (gmax > _gtol && iterg < _maxiterg) // Check for slip system resistance update tolerance
  {
    preSolveStress();
    solveStress();
    if (_err_tol)
      return;

    update_energies(); // Update plastic energy

    postSolveStress(); // Update _fp_old_inv = _fp_old

    gss_prev = _gss_tmp;

    update_slip_system_resistance(); // Update slip system resistance

    gmax = 0.0;
    for (unsigned i = 0; i < _nss; ++i)
    {
      gdiff = std::abs(gss_prev[i] - _gss_tmp[i]); // Calculate increment size

      if (gdiff > gmax)
        gmax = gdiff;
    }
    iterg++;
  }

  if (iterg == _maxiterg)
  {
#ifdef DEBUG
    mooseWarning("FiniteStrainCrystalPLasticity: Hardness Integration error gmax", gmax, "\n");
#endif
    _err_tol = true;
  }
}

void
FiniteStrainCrystalPlasticityMcAuliffe::postSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss[_qp] = _gss_tmp;
    _W0p[_qp] = _W0p_tmp;
    _acc_slip[_qp] = _accslip_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _gss[_qp] = _gss_tmp;
      _W0p[_qp] = _W0p_tmp;
      _acc_slip[_qp] = _accslip_tmp;
    }
    else
    {
      _gss_tmp_old = _gss_tmp;
      _W0p_tmp_old = _W0p_tmp;
      _accslip_tmp_old = _accslip_tmp;
    }
  }
}

// Update slip system resistance and elastic and plastic work
void
FiniteStrainCrystalPlasticityMcAuliffe::update_energies()
{
  RankTwoTensor cauchy_stress_undamaged, WpToTrace, invFe;
  Real detFe;
  Real c = _c[_qp];
  Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  Real kappa = lambda + (2.0 * mu / 3.0);

  if (_max_substep_iter == 1) //No substepping
  {
    _W0p_tmp = _W0p_old[_qp];
  }
  else
  {
    _W0p_tmp = _W0p_tmp_old;
  }

  // Update elastic and plastic work
  detFe = _fe.det();
  invFe = _fe.inverse();

  // _pk2[_qp] is the updated piola-kirchoff
  cauchy_stress_undamaged = _fe * _pk2_undamaged[_qp] * _fe.transpose()/detFe;
  WpToTrace = cauchy_stress_undamaged * _fe * ( _fp_inv.inverse() - _fp_old_inv.inverse() ) * _fp_inv * invFe * detFe;

  if ( detFe >= 0.0 ) {
    _W0e[_qp] = (kappa / 2.0) * ( 0.5 * (detFe*detFe - 1.0) - std::log(detFe) );
  }
  else {
    _W0e[_qp] = 0.0;
  }
  _W0p_tmp += WpToTrace.trace();
  //_dW0e_dstrain[_qp] = cauchy_stress_undamaged;
  _dW0p_dstrain[_qp].zero();
}

void
FiniteStrainCrystalPlasticityMcAuliffe::calcResidual( RankTwoTensor &resid )
{
  RankTwoTensor iden, ce, ce_pk2, eqv_slip_incr, pk2_new, pk2_new_dev, pk2_new_vol, invCe;
  Real trCe, detFe;
  Real c = _c[_qp];
  Real xfac = Utility::pow<2>(1.0-c) + _kdamage;

  // Isotropic elasticity is assumed
  Real lambda = _elasticity_tensor[_qp](0, 0, 1, 1);
  Real mu = _elasticity_tensor[_qp](0, 1, 0, 1);
  Real kappa = lambda + (2.0 * mu / 3.0);

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
  invCe = ce.inverse();
  trCe = ce.trace();
  detFe = _fe.det();

  pk2_new_vol = (kappa / 2.0) * ( detFe*detFe - 1.0 ) * invCe;
  pk2_new_dev = ( mu / std::pow(detFe,(2.0/3.0)) ) * ( iden - (trCe / 3.0) * invCe);
  if ( detFe >= 0.0 ) {
    pk2_new = pk2_new_vol * xfac;
  }
  else {
    pk2_new = pk2_new_vol;
  }
  pk2_new += pk2_new_dev;

  _dW0e_dstrain[_qp] = pk2_new_vol;

  // Used in StressDivergencePFFracTensors Jacobian
  if (c < 1.0)
    _dstress_dc[_qp] = -pk2_new_vol * (2.0 * (1.0 - c));
  else
    _dstress_dc[_qp].zero();

  _pk2_undamaged[_qp] = pk2_new_vol + pk2_new_dev;

  resid = _pk2_tmp - pk2_new;
}

// Calculate slip increment,dslipdtau. Override to modify.
void
FiniteStrainCrystalPlasticityMcAuliffe::getSlipIncrements()
{
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr(i) = _a0(i) * std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i)) *
                    copysign(1.0, _tau(i)) * _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt)
    {
      _slip_incr(i) = _slip_incr_tol * _dt * copysign(1.0, _tau(i));
      //_err_tol = true;
#ifdef DEBUG
      mooseWarning("Maximum allowable slip increment exceeded ", std::abs(_slip_incr(i)));
#endif
      return;
    }
  }

  for (unsigned int i = 0; i < _nss; ++i)
  {
    _dslipdtau(i) = _a0(i) / _xm(i) *
                    std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i) - 1.0) / _gss_tmp[i] *
                    _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt)
    {
      _dslipdtau(i) = 0.0;
    }
  }
}
