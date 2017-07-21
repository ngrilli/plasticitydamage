#include "plasticitydamageApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "CrystalPlasticitySlipRateDislo.h"
#include "CrystalPlasticitySlipResistanceDislo.h"
#include "CrystalPlasticityStateVarRateComponentDislo.h"
#include "NSMomentumInviscidFluxWithoutP.h"
#include "InertialForceCompressible.h"
#include "FiniteStrainCrystalPlasticityEOS.h"
#include "FiniteStrainCrystalPlasticityBulkViscosity.h"
#include "FiniteStrainCrystalPlasticityEOSBulkViscosity.h"
#include "FiniteStrainCrystalPlasticityEOSBulkViscosityFdot.h"
#include "FiniteStrainCrystalPlasticityEOSBulkViscositySigned.h"
#include "FiniteStrainCrystalPlasticityEOSquadC0.h"
#include "FiniteStrainCrystalPlasticityPhonon.h"
#include "FiniteStrainCrystalPlasticityDamage.h"
#include "FiniteStrainCrystalPlasticityMcAuliffe.h"
#include "FiniteStrainCrystalPlasticityW0p.h"
#include "FiniteStrainCrystalPlasticityDamagePrincipalStrains.h"
#include "LinearIsoElasticPFDamageEOS.h"
#include "LinearIsoElasticPFDamageEOSBulkViscosity.h"
#include "LinearIsoElasticPFAmor.h"
#include "LinearIsoElasticPFDamageNoc.h"
#include "LinearIsoElasticPFDamageMiehe.h"
#include "PFFracBulkRatePlastic.h"
#include "PFFracBulkRatePlasticFactor.h"
#include "PFFracBulkRateNob.h"
#include "PFFracBulkRateLaplacian.h"
#include "TimeDerivativeDamage.h"
#include "ReactionAbsolute.h"
#include "PFFracBulkRateNobLimited.h"
#include "PFFracBulkRateCAux.h"
#include "PFFracBulkRateMiehe.h"
#include "damageIC.h"
#include "damageICxyz.h"
#include "CrystalPlasticityPFDamageMiehe.h"
#include "CrackFrictionHeatSource.h"
#include "CrackPropagationHeatSource.h"
#include "PlasticHeatingSource.h"
#include "ComputeCrackFrictionHeatEnergy.h"
#include "ComputeCrackFrictionHeatEnergyDienes.h"
#include "ComputeExtraStressEOS.h"
#include "ThermalExpansionHeatSource.h"

template <>
InputParameters
validParams<plasticitydamageApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

plasticitydamageApp::plasticitydamageApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  plasticitydamageApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  plasticitydamageApp::associateSyntax(_syntax, _action_factory);
}

plasticitydamageApp::~plasticitydamageApp() {}

// External entry point for dynamic application loading
extern "C" void
plasticitydamageApp__registerApps()
{
  plasticitydamageApp::registerApps();
}
void
plasticitydamageApp::registerApps()
{
  registerApp(plasticitydamageApp);
}

// External entry point for dynamic object registration
extern "C" void
plasticitydamageApp__registerObjects(Factory & factory)
{
  plasticitydamageApp::registerObjects(factory);
}
void
plasticitydamageApp::registerObjects(Factory & factory)
{
  registerMaterial(FiniteStrainCrystalPlasticityEOS);
  registerMaterial(FiniteStrainCrystalPlasticityBulkViscosity);
  registerMaterial(FiniteStrainCrystalPlasticityEOSBulkViscosity);
  registerMaterial(FiniteStrainCrystalPlasticityEOSBulkViscosityFdot);
  registerMaterial(FiniteStrainCrystalPlasticityEOSBulkViscositySigned);
  registerMaterial(FiniteStrainCrystalPlasticityEOSquadC0);
  registerMaterial(FiniteStrainCrystalPlasticityPhonon);
  registerMaterial(FiniteStrainCrystalPlasticityDamage);
  registerMaterial(FiniteStrainCrystalPlasticityMcAuliffe);
  registerMaterial(FiniteStrainCrystalPlasticityW0p);
  registerMaterial(FiniteStrainCrystalPlasticityDamagePrincipalStrains);
  registerMaterial(LinearIsoElasticPFDamageEOS);
  registerMaterial(LinearIsoElasticPFDamageEOSBulkViscosity);
  registerMaterial(LinearIsoElasticPFAmor);
  registerMaterial(LinearIsoElasticPFDamageNoc);
  registerMaterial(LinearIsoElasticPFDamageMiehe);
  registerMaterial(CrystalPlasticityPFDamageMiehe);
  registerMaterial(ComputeCrackFrictionHeatEnergy);
  registerMaterial(ComputeCrackFrictionHeatEnergyDienes);
  registerMaterial(ComputeExtraStressEOS);

  registerKernel(NSMomentumInviscidFluxWithoutP);
  registerKernel(InertialForceCompressible);
  registerKernel(PFFracBulkRatePlastic);
  registerKernel(PFFracBulkRatePlasticFactor);
  registerKernel(PFFracBulkRateNob);
  registerKernel(PFFracBulkRateLaplacian);
  registerKernel(TimeDerivativeDamage);
  registerKernel(ReactionAbsolute);
  registerKernel(PFFracBulkRateNobLimited);
  registerKernel(PFFracBulkRateCAux);
  registerKernel(PFFracBulkRateMiehe);
  registerKernel(CrackFrictionHeatSource);
  registerKernel(CrackPropagationHeatSource);
  registerKernel(PlasticHeatingSource);
  registerKernel(ThermalExpansionHeatSource);

  registerUserObject(CrystalPlasticitySlipRateDislo);
  registerUserObject(CrystalPlasticitySlipResistanceDislo);
  registerUserObject(CrystalPlasticityStateVarRateComponentDislo);

  registerInitialCondition(damageIC);
  registerInitialCondition(damageICxyz);
}

// External entry point for dynamic syntax association
extern "C" void
plasticitydamageApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  plasticitydamageApp::associateSyntax(syntax, action_factory);
}
void
plasticitydamageApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
