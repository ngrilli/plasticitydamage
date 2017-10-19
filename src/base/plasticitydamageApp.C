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
#include "FiniteStrainCrystalPlasticityPhonon.h"
#include "FiniteStrainCrystalPlasticityDamage.h"
#include "FiniteStrainCrystalPlasticityMcAuliffe.h"
#include "FiniteStrainCrystalPlasticityW0p.h"
#include "FiniteStrainCrystalPlasticityDamagePrincipalStrains.h"
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
#include "CrackPropagationHeatSourceNoDevel.h"
#include "PlasticHeatingSource.h"
#include "ComputeArrheniusMassFractionRateLimit.h"
#include "ComputeArrheniusMassFractionRateLimit3species.h"
#include "ArrheniusHeatEnergy.h"
#include "ArrheniusMassFraction.h"
#include "ArrheniusMassFractionRateLimit.h"
#include "ArrheniusMassFractionRateLimit3species.h"
#include "ArrheniusHeatEnergyRateLimit.h"
#include "PlasticHeatingSourceMiehe2016.h"
#include "ComputeCrackFrictionHeatEnergy.h"
#include "ComputeCrackFrictionHeatEnergyDienes.h"
#include "ComputeCrackFrictionHeatEnergyDienesFiniteStrain.h"
#include "ThermalExpansionHeatSource.h"
#include "PFThermalConductivityCompressive.h"
#include "HeatConductionMaterial3species.h"

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
  registerMaterial(FiniteStrainCrystalPlasticityPhonon);
  registerMaterial(FiniteStrainCrystalPlasticityDamage);
  registerMaterial(FiniteStrainCrystalPlasticityMcAuliffe);
  registerMaterial(FiniteStrainCrystalPlasticityW0p);
  registerMaterial(FiniteStrainCrystalPlasticityDamagePrincipalStrains);
  registerMaterial(LinearIsoElasticPFAmor);
  registerMaterial(LinearIsoElasticPFDamageNoc);
  registerMaterial(LinearIsoElasticPFDamageMiehe);
  registerMaterial(CrystalPlasticityPFDamageMiehe);
  registerMaterial(ComputeCrackFrictionHeatEnergy);
  registerMaterial(ComputeCrackFrictionHeatEnergyDienes);
  registerMaterial(ComputeCrackFrictionHeatEnergyDienesFiniteStrain);
  registerMaterial(PFThermalConductivityCompressive);
  registerMaterial(ComputeArrheniusMassFractionRateLimit);
  registerMaterial(ComputeArrheniusMassFractionRateLimit3species);
  registerMaterial(HeatConductionMaterial3species);

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
  registerKernel(CrackPropagationHeatSourceNoDevel);
  registerKernel(PlasticHeatingSource);
  registerKernel(ArrheniusHeatEnergy);
  registerKernel(ArrheniusMassFraction);
  registerKernel(ArrheniusMassFractionRateLimit);
  registerKernel(ArrheniusMassFractionRateLimit3species);
  registerKernel(ArrheniusHeatEnergyRateLimit);
  registerKernel(PlasticHeatingSourceMiehe2016);
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
