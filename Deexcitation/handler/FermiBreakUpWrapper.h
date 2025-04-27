#pragma once

#include <memory>

#include "G4VFermiBreakUp.hh"
#include "FermiBreakUp/FermiBreakUp.h"

class FermiBreakUpWrapper : public G4VFermiBreakUp {
 public:
  FermiBreakUpWrapper();
  FermiBreakUpWrapper(fbu::FermiBreakUp&& model);

  void Initialise() override;

  void BreakFragment(G4FragmentVector* fragmentsPtr, G4Fragment* fragment) override;

  static G4bool IsFermiPossible(G4int Z, G4int A, G4double excitationEnergy);

  G4bool IsApplicable(G4int Z, G4int A, G4double excitationEnergy) const override;

 private:
  fbu::FermiBreakUp fermiModel_;
};
