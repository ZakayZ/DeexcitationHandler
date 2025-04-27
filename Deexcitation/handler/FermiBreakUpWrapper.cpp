//
// Created by Artem Novikov on 20.05.2023.
//

#include <G4BaryonConstructor.hh>
#include <memory>

#include "FermiBreakUp/fragment_pool/FragmentsStorage.h"
#include "FermiBreakUp/fragment_pool/FragmentPool.h"
#include "FermiBreakUp/nuclei_properties/NucleiProperties.h"
#include "FermiBreakUp/util/DataTypes.h"
#include "FermiBreakUp/util/Cache.h"
#include "G4Fragment.hh"
#include "FermiBreakUpWrapper.h"

constexpr std::int32_t MAX_A = 19;
constexpr std::int32_t MAX_Z = 9;

FermiBreakUpWrapper::FermiBreakUpWrapper()
  : fermiModel_(fbu::FermiBreakUp(std::make_unique<fbu::LFUCache<fbu::NucleiData, fbu::FragmentSplits>>(MAX_A * MAX_A / 2)))
{
}

FermiBreakUpWrapper::FermiBreakUpWrapper(fbu::FermiBreakUp&& model)
  : fermiModel_(std::move(model))
{
}

void FermiBreakUpWrapper::BreakFragment(G4FragmentVector* fragmentsPtr, G4Fragment* fragment) {
  auto results = fermiModel_.BreakItUp(fbu::Particle(fbu::AtomicMass(fragment->GetA_asInt()),
                                                 fbu::ChargeNumber(fragment->GetZ_asInt()),
                                                 fragment->GetMomentum()));

  for (const auto& particle : results) {
    fragmentsPtr->push_back(new G4Fragment(G4int(particle.GetAtomicMass()),
                                           G4int(particle.GetChargeNumber()),
                                           G4LorentzVector(particle.GetMomentum())));
    fragmentsPtr->back()->SetCreatorModelID(24300);
    fragmentsPtr->back()->SetCreationTime(fragment->GetCreationTime());
  }
}

void FermiBreakUpWrapper::Initialise() {
  // order matters, FragmentPool uses NucleiProperties
  fbu::NucleiProperties::Reset();
  if (G4NucleiProperties::GetNuclearMass(2, 0) <= 0.) {
    G4BaryonConstructor pCBar;
    pCBar.ConstructParticle();
  }
  for (auto a = 1; a < MAX_A; ++a) {
    for (auto z = 0; z <= a; ++z) {
      const auto atomicMass = fbu::AtomicMass(a);
      const auto chargeNumber = fbu::ChargeNumber(z);

      const auto mass = G4NucleiProperties::GetNuclearMass(a, z);
      if (mass > 0.) {
        fbu::NucleiProperties::Instance().AddNuclei(atomicMass, chargeNumber, mass, G4NucleiProperties::IsInStableTable(a, z));
      }
    }
  }

  fbu::FragmentPool::Reset();
}

G4bool FermiBreakUpWrapper::IsFermiPossible(G4int Z, G4int A, [[maybe_unused]] G4double excitationEnergy) {
  return Z < MAX_Z && A < MAX_A;  // && excitationEnergy > -10 * CLHEP::keV;
}

G4bool FermiBreakUpWrapper::IsApplicable(G4int Z, G4int A, G4double excitationEnergy) const {
  return IsFermiPossible(Z, A, excitationEnergy);
}
