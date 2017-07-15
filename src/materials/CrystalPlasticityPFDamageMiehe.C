/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CrystalPlasticityPFDamageMiehe.h"
#include "petscblaslapack.h"
#include "libmesh/utility.h"
#include "MathUtils.h"

template<>
InputParameters validParams<CrystalPlasticityPFDamageMiehe>()
{
  InputParameters params = validParams<FiniteStrainCrystalPlasticity>();
  params.addClassDescription("Crystal Plasticity base class: FCC system with power law flow rule implemented. Damage. Miehe 2016");
  params.addRequiredCoupledVar("b", "Laplacian of the order parameter for damage");
  params.addRequiredCoupledVar("initialc", "initial damage");
  params.addRequiredParam<Real>("l", "Interface width");
  params.addRequiredParam<Real>("visco", "Viscosity parameter");
  params.addRequiredParam<Real>("Wc","Threshold fracture energy");
  params.addRequiredParam<MaterialPropertyName>("gc_prop_var","Material property name with gc value");
  params.addParam<Real>("kdamage",1e-6,"Stiffness of damaged matrix");
  params.addRequiredParam<Real>("n_Murnaghan", "exponent in Murnaghan EOS");
  params.addRequiredParam<Real>("bulk_modulus_ref", "reference bulk modulus");
  params.addRequiredParam<Real>("C0", "Von Neuman coefficient");
  params.addRequiredParam<Real>("C1", "Landshoff coefficient");

  return params;
}

CrystalPlasticityPFDamageMiehe::CrystalPlasticityPFDamageMiehe(const InputParameters & parameters) :
    FiniteStrainCrystalPlasticity(parameters),
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
    _x_bracket(declareProperty<Real>("x_bracket")),
    _x_bracket_old(declarePropertyOld<Real>("x_bracket")),
    _l(getParam<Real>("l")),
    _visco(getParam<Real>("visco")),
    _Wc(getParam<Real>("Wc")),
    _n_Murnaghan(getParam<Real>("n_Murnaghan")),
    _Bulk_Modulus_Ref(getParam<Real>("bulk_modulus_ref")),
    _C0(getParam<Real>("C0")),
    _C1(getParam<Real>("C1")),
    _W0e(declareProperty<Real>("W0e")), // elastic energy
    _W0e_old(declarePropertyOld<Real>("W0e")), // elastic energy of previous increment
    _W0p(declareProperty<Real>("W0p")), // plastic energy
    _W0p_old(declarePropertyOld<Real>("W0p")), // plastic energy of previous increment
    //_dstress_dc(declarePropertyDerivative<RankTwoTensor>(_base_name + "stress", getVar("c", 0)->name())),
    _dW0e_dstrain(declareProperty<RankTwoTensor>("dW0e_dstrain")),
    _dW0p_dstrain(declareProperty<RankTwoTensor>("dW0p_dstrain")),
    _fe_old(declareProperty<RankTwoTensor>("fe_old")),
    _cauchy_out(declareProperty<RankTwoTensor>("cauchy_out"))
{
}

/**
 * Solves stress residual equation using NR.
 * Updates slip system resistances iteratively.
 */
void
CrystalPlasticityPFDamageMiehe::computeQpStress()
{
  unsigned int substep_iter = 1; // Depth of substepping; Limited to maximum substep iteration
  unsigned int num_substep = 1;  // Calculated from substep_iter as 2^substep_iter
  Real dt_original = _dt;        // Stores original _dt; Reset at the end of solve
  _first_substep = true;         // Initialize variables at substep_iter = 1

  if (_max_substep_iter > 1)
  {
    _dfgrd_tmp_old = _deformation_gradient_old[_qp];
    if (_dfgrd_tmp_old.det() == 0)
      _dfgrd_tmp_old.addIa(1.0);

    _delta_dfgrd = _deformation_gradient[_qp] - _dfgrd_tmp_old;
    _err_tol = true; // Indicator to continue substepping
  }

  // Substepping loop
  while (_err_tol && _max_substep_iter > 1)
  {
    _dt = dt_original / num_substep;

    for (unsigned int istep = 0; istep < num_substep; ++istep)
    {
      _first_step_iter = false;
      if (istep == 0)
        _first_step_iter = true;

      _last_step_iter = false;
      if (istep == num_substep - 1)
        _last_step_iter = true;

      _dfgrd_scale_factor = (static_cast<Real>(istep) + 1) / num_substep;
      _dfgrd_scale_factor_old = (static_cast<Real>(istep)) / num_substep;

      preSolveQp();
      solveQp();

      if (_err_tol)
      {
        substep_iter++;
        num_substep *= 2;
        break;
      }
    }

    _first_substep = false; // Prevents reinitialization
    _dt = dt_original;      // Resets dt

#ifdef DEBUG
    if (substep_iter > _max_substep_iter)
      mooseWarning("FiniteStrainCrystalPlasticity: Failure with substepping");
#endif

    if (!_err_tol || substep_iter > _max_substep_iter)
      postSolveQp(); // Evaluate variables after successful solve or indicate failure
  }

  // No substepping
  if (_max_substep_iter == 1)
  {
    preSolveQp();
    solveQp();
    postSolveQp();
  }
}

void
CrystalPlasticityPFDamageMiehe::preSolveQp()
{
  // Initialize variable
  if (_first_substep)
  {
    _Jacobian_mult[_qp].zero(); // Initializes jacobian for preconditioner
    calc_schmid_tensor();
  }

  if (_max_substep_iter == 1)
    {
      _dfgrd_tmp = _deformation_gradient[_qp]; // Without substepping
      _dfgrd_tmp_old_substep = _deformation_gradient_old[_qp];
    }
  else
    {
      _dfgrd_tmp = _dfgrd_scale_factor * _delta_dfgrd + _dfgrd_tmp_old;
      _dfgrd_tmp_old_substep = _dfgrd_scale_factor_old * _delta_dfgrd + _dfgrd_tmp_old;
    }

  _err_tol = false;
}

void
CrystalPlasticityPFDamageMiehe::preSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss_tmp = _gss_old[_qp];
    _W0e_tmp = _W0e_old[_qp];
    _W0p_tmp = _W0p_old[_qp];
    _accslip_tmp_old = _acc_slip_old[_qp];
    _damage_tmp_old = _damage_old[_qp];
    _b_tmp_old = _b_old[_qp];
  }
  else
  {
    if (_first_step_iter)
    {
      _gss_tmp = _gss_tmp_old = _gss_old[_qp];
      _W0e_tmp = _W0e_tmp_old = _W0e_old[_qp];
      _W0p_tmp = _W0p_tmp_old = _W0p_old[_qp];
      _accslip_tmp = _accslip_tmp_old = _acc_slip_old[_qp];
      _damage_tmp = _damage_tmp_old = _damage_old[_qp];
      _b_tmp = _b_tmp_old = _b_old[_qp];
    }
    else
    {
      _gss_tmp = _gss_tmp_old;
      _W0e_tmp = _W0e_tmp_old;
      _W0p_tmp = _W0p_tmp_old;
      _damage_tmp = _damage_tmp_old;
      _b_tmp = _b_tmp_old;
    }
  }
}

void
CrystalPlasticityPFDamageMiehe::solveStatevar()
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

    update_energies();

    update_damage();

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
CrystalPlasticityPFDamageMiehe::postSolveStatevar()
{
  if (_max_substep_iter == 1)//No substepping
  {
    _gss[_qp] = _gss_tmp;
    _W0e[_qp] = _W0e_tmp;
    _W0p[_qp] = _W0p_tmp;
    _acc_slip[_qp] = _accslip_tmp;
    _damage[_qp] = _damage_tmp;
  }
  else
  {
    if (_last_step_iter)
    {
      _gss[_qp] = _gss_tmp;
      _W0e[_qp] = _W0e_tmp;
      _W0p[_qp] = _W0p_tmp;
      _acc_slip[_qp] = _accslip_tmp;
      _damage[_qp] = _damage_tmp;
    }
    else
    {
      _gss_tmp_old = _gss_tmp;
      _W0e_tmp_old = _W0e_tmp;
      _W0p_tmp_old = _W0p_tmp;
      _accslip_tmp_old = _accslip_tmp;
      _damage_tmp_old = _damage_tmp;
      _b_tmp_old = _b_tmp;
    }
  }
}

// Update slip system resistance and elastic and plastic work
void
CrystalPlasticityPFDamageMiehe::update_energies()
{
  RankTwoTensor cauchy_stress, WeToTrace, WpToTrace, invFe;
  Real detFe, detFe_old;
  Real c_tmp_old = _initialc_old[_qp] + _damage_tmp_old;

  // Update elastic and plastic work
  detFe = _fe.det();
  _fe_tmp_old = _dfgrd_tmp_old_substep * _fp_old_inv;
  _fe_old[_qp] = _fe_tmp_old;
  detFe_old = _fe_tmp_old.det();
  invFe = _fe.inverse();

  // _pk2[_qp] is the updated piola-kirchoff
  cauchy_stress = _fe * _pk2[_qp] * _fe.transpose()/detFe;
  _cauchy_out[_qp] = cauchy_stress;
  WeToTrace = cauchy_stress * ( _fe - _fe_tmp_old ) * invFe * detFe;
  WpToTrace = cauchy_stress * _fe * ( _fp_inv.inverse() - _fp_old_inv.inverse() ) * _fp_inv * invFe * detFe;
  _W0e_tmp += WeToTrace.trace();
  _W0p_tmp += WpToTrace.trace();
  _dW0e_dstrain[_qp] = cauchy_stress;
  _dW0p_dstrain[_qp] = cauchy_stress;
  //_dstress_dc[_qp] = -cauchy_stress * (2.0 * (1.0 - c_tmp_old));

}

// Update damage
void
CrystalPlasticityPFDamageMiehe::update_damage()
{
  Real c_tmp_old = _initialc_old[_qp] + _damage_tmp_old;
  const Real gc = _gc_prop[_qp];
  Real signx;
  Real xfac = _kdamage;
  if (c_tmp_old < 1.0)
    xfac += Utility::pow<2>(1.0 - c_tmp_old);

  _x_bracket[_qp] = _l*_b_tmp_old + 2.0*(1.0-c_tmp_old)*(_W0e_tmp+_W0p_tmp-_Wc)/gc - c_tmp_old/_l;
  signx = MathUtils::sign(_x_bracket[_qp]);
  _x_bracket[_qp] = ((signx + 1.0) / 2.0) * _x_bracket[_qp] / _visco;

  _damage_tmp = _damage_tmp_old + _x_bracket[_qp] * _dt;

  if ( ( _damage_tmp + _initialc_old[_qp] ) > 1.0 ) {
    _damage_tmp = 1.0 - _initialc_old[_qp];
  }
  _c[_qp] = _damage_tmp + _initialc_old[_qp];

}

void
CrystalPlasticityPFDamageMiehe::calcResidual( RankTwoTensor & resid )
{
  RankTwoTensor iden, ce, ee, ce_pk2, eqv_slip_incr, pk2_new;
  Real trD;
  Real c = _initialc_old[_qp] + _damage_tmp_old;
  Real xfac = Utility::pow<2>(1.0-c) + _kdamage;

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
  ee = ce - iden;
  ee *= 0.5;

  pk2_new = _elasticity_tensor[_qp] * ee;

  pk2_new.addIa(-1.0/3.0 * pk2_new.trace());
  pk2_new.addIa( (_Bulk_Modulus_Ref/_n_Murnaghan) * ( 1.0 - std::pow( 1.0/_fe.det() , _n_Murnaghan ) ) );

  pk2_new = xfac * pk2_new;

  trD = ( _dfgrd_tmp.det() - _dfgrd_tmp_old_substep.det() ) / _dt;
  trD /= _dfgrd_tmp_old_substep.det();

  pk2_new.addIa( _C0 * trD * trD );
  pk2_new.addIa( _C1 * abs(trD) );

  resid = _pk2_tmp - pk2_new;
}

// Calculate slip increment,dslipdtau. Override to modify.
void
CrystalPlasticityPFDamageMiehe::getSlipIncrements()
{
  for (unsigned int i = 0; i < _nss; ++i)
  {
    _slip_incr(i) = _a0(i) * std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i)) * copysign(1.0, _tau(i)) * _dt;
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
    _dslipdtau(i) = _a0(i) / _xm(i) * std::pow(std::abs(_tau(i) / _gss_tmp[i]), 1.0 / _xm(i) - 1.0) / _gss_tmp[i] * _dt;
    if (std::abs(_slip_incr(i)) > _slip_incr_tol * _dt)
    {
      _dslipdtau(i) = 0.0;
    }
  }
}
