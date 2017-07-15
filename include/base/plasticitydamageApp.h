#ifndef PLASTICITYDAMAGEAPP_H
#define PLASTICITYDAMAGEAPP_H

#include "MooseApp.h"

class plasticitydamageApp;

template <>
InputParameters validParams<plasticitydamageApp>();

class plasticitydamageApp : public MooseApp
{
public:
  plasticitydamageApp(InputParameters parameters);
  virtual ~plasticitydamageApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* PLASTICITYDAMAGEAPP_H */
