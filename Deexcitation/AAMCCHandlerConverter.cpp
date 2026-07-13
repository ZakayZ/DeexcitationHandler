// NOLINTNEXTLINE(misc-include-cleaner)
#include "Deexcitation/AAMCCHandlerConverter.h"

#include "Deexcitation/aamcc/include/DeexcitationHandler.hh"

#include <G4Fragment.hh>
#include <G4ReactionProduct.hh>

namespace {

  G4Fragment ColaToG4(const cola::Particle& particle) {
    const auto [atomic_mass, charge] = particle.GetAZ();

    return G4Fragment(
        static_cast<G4int>(atomic_mass), static_cast<G4int>(charge),
        G4LorentzVector(particle.momentum.x, particle.momentum.y, particle.momentum.z, particle.momentum.e));
  }

  cola::Particle G4ToCola(const G4ReactionProduct& fragment) {
    return cola::Particle{
        .position = cola::LorentzVector{.e = 0., .x = 0., .y = 0., .z = 0.},
        .momentum =
            cola::LorentzVector{
                .e = fragment.GetTotalEnergy(),
                .x = fragment.GetMomentum().x(),
                .y = fragment.GetMomentum().y(),
                .z = fragment.GetMomentum().z(),
            },
        .pdg_code = cola::AZToPdg({static_cast<uint16_t>(fragment.GetDefinition()->GetAtomicMass()),
                                   static_cast<uint16_t>(fragment.GetDefinition()->GetAtomicNumber())}),
        .p_class = cola::ParticleClass::kProduced,
    };
  }

}  // namespace

namespace cola {

  AAMCCHandlerConverter::AAMCCHandlerConverter(std::unique_ptr<DeexcitationHandler>&& model)
      : model_(std::move(model)) {}

  std::unique_ptr<EventData> AAMCCHandlerConverter::operator()(std::unique_ptr<EventData>&& data) {
    EventParticles results;
    for (const auto& particle : data->particles) {
      const auto particle_class = particle.p_class;
      if (particle_class == ParticleClass::kSpectatorA || particle_class == ParticleClass::kSpectatorB) {
        G4ReactionProductVector* model_result = model_->AAMCCBreakItUp(ColaToG4(particle));
        if (model_result != nullptr) {
          for (const auto* fragment : *model_result) {
            if (fragment != nullptr) {
              results.emplace_back(G4ToCola(*fragment));
              results.back().p_class = particle_class;
            }
          }
          for (auto* fragment : *model_result) {
            delete fragment;
          }
          delete model_result;
        }
      } else {
        results.push_back(particle);
      }
    }

    data->particles = std::move(results);
    return std::move(data);
  }

}  // namespace cola
