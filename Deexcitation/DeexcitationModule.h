#ifndef DEEXCITATION_DEEXCITATIONMODULE_H_
#define DEEXCITATION_DEEXCITATIONMODULE_H_

#include "Deexcitation/AAMCCHandlerFactory.h"
#include "Deexcitation/G4HandlerFactory.h"

#include <COLA.hh>

namespace cola {

  using DeexcitationModule = GenericModule<G4HandlerFactory, AAMCCHandlerFactory>;

}  // namespace cola

#endif  // DEEXCITATION_DEEXCITATIONMODULE_H_
