//
// Created by Artem Novikov on 17.05.2023.
//

#include <string_view>

#include <CLHEP/Units/PhysicalConstants.h>

#include <G4BosonConstructor.hh>
#include <G4LeptonConstructor.hh>
#include <G4MesonConstructor.hh>
#include <G4BaryonConstructor.hh>
#include <G4IonConstructor.hh>
#include <G4ProcessManager.hh>
#include <G4StateManager.hh>
#include <G4RunManager.hh>

#include <G4LorentzVector.hh>
#include <G4NistManager.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleTypes.hh>
#include <G4Ions.hh>
#include <G4Electron.hh>

#include <G4Evaporation.hh>
#include <G4PhotonEvaporation.hh>
#include <G4StatMF.hh>

#include "FermiBreakUpWrapper.h"
#include "ExcitationHandler.h"

namespace {
  constexpr size_t EvaporationIterationThreshold = 1e3;

  static const std::string ErrorNoModel = "no model was applied, check conditions";

  void CleanUp(G4FragmentVector& v, std::queue<G4Fragment*>& q1, std::queue<G4Fragment*>& q2) {
    for (auto ptr : v) {
      delete ptr;
    }
    while (!q1.empty()) {
      delete q1.front();
      q1.pop();
    }
    while (!q2.empty()) {
      delete q2.front();
      q2.pop();
    }
  }

  constexpr G4int HashParticle(G4int A, G4int Z) { return A * 1000 + Z; }

  G4ParticleDefinition* SpecialParticleDefinition(const G4Fragment& fragment) {
    switch (HashParticle(fragment.GetA_asInt(), fragment.GetZ_asInt())) {
      case HashParticle(0, 0): {
        return G4Gamma::GammaDefinition();
      }

      case HashParticle(0, -1): {
        return G4Electron::ElectronDefinition();
      }

      case HashParticle(-1, 1): {
        return G4PionPlus::PionPlus();
      }

      case HashParticle(-1, -1): {
        return G4PionMinus::PionMinus();
      }

      case HashParticle(-1, 0): {
        return G4PionZero::PionZero();
      }

      case HashParticle(1, 0): {
        return G4Neutron::NeutronDefinition();
      }

      case HashParticle(1, 1): {
        return G4Proton::ProtonDefinition();
      }

      case HashParticle(2, 1): {
        return G4Deuteron::DeuteronDefinition();
      }

      case HashParticle(3, 1): {
        return G4Triton::TritonDefinition();
      }

      case HashParticle(3, 2): {
        return G4He3::He3Definition();
      }

      case HashParticle(4, 2): {
        return G4Alpha::AlphaDefinition();
      }
    }

    return nullptr;
  }

  void EvaporationError(const G4Fragment& fragment, const G4Fragment& currentFragment, size_t iter) {
    G4ExceptionDescription ed;
    ed << "Infinite loop in the de-excitation module: " << iter
      << " iterations \n"
      << "      Initial fragment: \n" << fragment
      << "\n      Current fragment: \n" << currentFragment;
    G4Exception("ExcitationHandler::BreakItUp", "", FatalException,
                ed, "Stop execution");
  }
} // namespace


ExcitationHandler::ExcitationHandler()
  : multiFragmentationModel_(DefaultMultiFragmentation())
  , fermiBreakUpModel_(DefaultFermiBreakUp())
  , photonEvaporationModel_(DefaultPhotonEvaporation())
  , evaporationModel_(DefaultEvaporation())
  , multiFragmentationCondition_(DefaultMultiFragmentationCondition())
  , fermiCondition_(DefaultFermiBreakUpCondition())
  , photonEvaporationCondition_(DefaultPhotonEvaporationCondition())
  , evaporationCondition_(DefaultEvaporationCondition())
{
  evaporationModel_->SetFermiBreakUp(fermiBreakUpModel_.get());
  evaporationModel_->SetPhotonEvaporation(photonEvaporationModel_.get());

  G4BosonConstructor pCBos;
  pCBos.ConstructParticle();

  G4LeptonConstructor pCLept;
  pCLept.ConstructParticle();

  G4MesonConstructor pCMes;
  pCMes.ConstructParticle();

  G4BaryonConstructor pCBar;
  pCBar.ConstructParticle();

  G4IonConstructor pCIon;
  pCIon.ConstructParticle();

  G4GenericIon* gion = G4GenericIon::GenericIon();
  auto manager = new G4ProcessManager(gion);
  manager->SetVerboseLevel(0);
  gion->SetProcessManager(manager);

  G4StateManager::GetStateManager()->SetNewState(G4State_Init); // To let create ions
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ionTable = particleTable->GetIonTable();
  particleTable->SetReadiness();
  ionTable->CreateAllIon();
  ionTable->CreateAllIsomer();
}

ExcitationHandler::~ExcitationHandler() {
  photonEvaporationModel_.release();  // otherwise, SegFault in evaporation destructor
}

std::vector<G4ReactionProduct> ExcitationHandler::BreakItUp(const G4Fragment& fragment) {
  auto nist = G4NistManager::Instance();
  G4FragmentVector results;
  std::queue<G4Fragment*> evaporationQueue;
  std::queue<G4Fragment*> photonEvaporationQueue;

  // In case A <= 1 the fragment will not perform any nucleon emission
  auto initialFragmentPtr = std::make_unique<G4Fragment>(fragment);
  if (IsStable(fragment, nist)) {
    results.push_back(initialFragmentPtr.release());
  } else {
    if (multiFragmentationCondition_(fragment)) {
      ApplyMultiFragmentation(std::move(initialFragmentPtr), results, evaporationQueue);
    } else {
      evaporationQueue.push(initialFragmentPtr.release());
    }

    for (size_t iterationCount = 0; !evaporationQueue.empty(); ++iterationCount) {
      auto fragmentPtr = std::unique_ptr<G4Fragment>(evaporationQueue.front());
      evaporationQueue.pop();

      // infinite loopQEQ
      if (iterationCount == EvaporationIterationThreshold) {
        CleanUp(results, evaporationQueue, photonEvaporationQueue);

        EvaporationError(fragment, *fragmentPtr, iterationCount);
        return {};
        // process is dead
      }

      // FermiBreakUp part
      if (fermiCondition_(*fragmentPtr)) {
        ApplyFermiBreakUp(std::move(fragmentPtr), results, photonEvaporationQueue);
        continue;
      }

      // Evaporation part
      if (evaporationCondition_(*fragmentPtr)) {
        ApplyEvaporation(std::move(fragmentPtr), results, evaporationQueue);
        continue;
      }

      CleanUp(results, evaporationQueue, photonEvaporationQueue);
      throw std::runtime_error(ErrorNoModel);
    }

    // Photon EvapPation part
    while (!photonEvaporationQueue.empty()) {
      auto fragmentPtr = std::unique_ptr<G4Fragment>(photonEvaporationQueue.front());
      photonEvaporationQueue.pop();

      if (photonEvaporationCondition_(*fragmentPtr)) {
        ApplyPhotonEvaporation(std::move(fragmentPtr), results);
        continue;
      }

      CleanUp(results, evaporationQueue, photonEvaporationQueue);
      throw std::runtime_error(ErrorNoModel);
    }
  }

  auto reactionProducts = ConvertResults(results);

  CleanUp(results, evaporationQueue, photonEvaporationQueue);

  return reactionProducts;
}

std::unique_ptr<G4VMultiFragmentation> ExcitationHandler::DefaultMultiFragmentation() {
  return std::make_unique<G4StatMF>();
}

std::unique_ptr<G4VFermiBreakUp> ExcitationHandler::DefaultFermiBreakUp() {
  auto model = std::make_unique<FermiBreakUpWrapper>();
  model->Initialise();
  return model;
}

std::unique_ptr<G4VEvaporation> ExcitationHandler::DefaultEvaporation() {
  auto evaporation = std::make_unique<G4Evaporation>();
  return evaporation;
}

std::unique_ptr<G4VEvaporationChannel> ExcitationHandler::DefaultPhotonEvaporation() {
  return std::make_unique<G4PhotonEvaporation>();
}

ExcitationHandler::Condition ExcitationHandler::DefaultMultiFragmentationCondition() {
  return [](const G4Fragment& fragment) -> bool {
    constexpr G4int maxAtomicMass = 19;
    constexpr G4int maxCharge = 9;
    constexpr G4double lowerBoundTransitionMF = 3 * CLHEP::MeV;
    constexpr G4double upperBoundTransitionMF = 5 * CLHEP::MeV;

    const auto atomicMass = fragment.GetA_asInt();
    const auto charge = fragment.GetZ_asInt();

    if (atomicMass < maxAtomicMass && charge < maxCharge) { return false; }

    const auto exitationEnergy = fragment.GetExcitationEnergy();

    const auto scale = 1. / (2. * (upperBoundTransitionMF - lowerBoundTransitionMF));
    const auto energyOffset = (upperBoundTransitionMF + lowerBoundTransitionMF) / 2.;
    const auto transitionProb = 0.5 * std::tanh((exitationEnergy / atomicMass - energyOffset) / scale) + 0.5;

    const auto random = G4RandFlat::shoot();

    if (exitationEnergy < lowerBoundTransitionMF * atomicMass) { return false; }

    if (random < transitionProb && exitationEnergy < upperBoundTransitionMF * atomicMass) { return true; }

    if (random > transitionProb && exitationEnergy < upperBoundTransitionMF * atomicMass) { return false; }

    if (exitationEnergy > upperBoundTransitionMF * atomicMass) { return true; }

    return false;
  };
}

ExcitationHandler::Condition ExcitationHandler::DefaultFermiBreakUpCondition() {
  return [](const G4Fragment& fragment) {
    return FermiBreakUpWrapper::IsFermiPossible(fragment.GetZ_asInt(),
                                              fragment.GetA_asInt(),
                                              fragment.GetExcitationEnergy());
  };
}

ExcitationHandler::Condition ExcitationHandler::DefaultEvaporationCondition() {
  return [](const G4Fragment&) { return true; };
}

ExcitationHandler::Condition ExcitationHandler::DefaultPhotonEvaporationCondition() {
  return [](const G4Fragment&) { return true; };
}

bool ExcitationHandler::IsGroundState(const G4Fragment& fragment) const {
  return fragment.GetExcitationEnergy() < stableThreshold_;
}

bool ExcitationHandler::IsStable(const G4Fragment& fragment, const G4NistManager* nist) const {
  return fragment.GetA_asInt() <= 1
         || (IsGroundState(fragment)
             && nist->GetIsotopeAbundance(fragment.GetZ_asInt(), fragment.GetA_asInt()) > 0);
}

void ExcitationHandler::ApplyMultiFragmentation(std::unique_ptr<G4Fragment> fragment,
                                                G4FragmentVector& results,
                                                std::queue<G4Fragment*>& nextStage) {
  auto fragments = std::unique_ptr<G4FragmentVector>(multiFragmentationModel_->BreakItUp(*fragment));
  if (fragments == nullptr || fragments->size() <= 1) {
    nextStage.push(fragment.release());
    return;
  }

  GroupFragments(*fragments, results, nextStage);
}

void ExcitationHandler::ApplyFermiBreakUp(std::unique_ptr<G4Fragment> fragment,
                                          G4FragmentVector& results,
                                          std::queue<G4Fragment*>& nextStage) {
  G4FragmentVector fragments;
  fermiBreakUpModel_->BreakFragment(&fragments, fragment.get());
  // auto fragments = std::unique_ptr<G4FragmentVector>(fermiBreakUpModel_->BreakItUp(fragment.get()))

  if (fragments.size() == 1) {
    nextStage.emplace(fragment.release());
    return;
  }

  GroupFragments(fragments, results, nextStage);
}

void ExcitationHandler::ApplyEvaporation(std::unique_ptr<G4Fragment> fragment,
                                         G4FragmentVector& results,
                                         std::queue<G4Fragment*>& nextStage) {
  G4FragmentVector fragments;
  evaporationModel_->BreakFragment(&fragments, fragment.get());
  // auto fragments = std::unique_ptr<G4FragmentVector>(evaporationModel_->BreakItUp(fragment.get()))

  // because evaporation adjusts it
  auto fragmentPtr = fragment.release();
  if (fragments.empty() || fragments.back() != fragmentPtr) {
    fragments.emplace_back(fragmentPtr);
  }

  if (fragments.size() == 1) {
    results.emplace_back(fragmentPtr);
    return;
  }

  GroupFragments(fragments, results, nextStage);
}

void ExcitationHandler::ApplyPhotonEvaporation(std::unique_ptr<G4Fragment> fragment, G4FragmentVector& results) {
  // photon de-excitation only for hot fragments
  if (!IsGroundState(*fragment)) {
    G4FragmentVector fragments;

    photonEvaporationModel_->BreakUpChain(&fragments, fragment.get());

    for (auto fragmentPtr : fragments) {
      results.emplace_back(fragmentPtr);
    }
  }

  // primary fragment is kept
  results.emplace_back(fragment.release());
}

void ExcitationHandler::GroupFragments(const G4FragmentVector& fragments,
                                       G4FragmentVector& results,
                                       std::queue<G4Fragment*>& nextStage) {
  auto nist = G4NistManager::Instance();

  // fragment pointers is moved to unique and will be deleted later
  for (auto fragmentPtr : fragments) {
    // gamma, p, n or stable nuclei
    if (IsStable(*fragmentPtr, nist)) {
      results.emplace_back(fragmentPtr);
    } else {
      nextStage.emplace(fragmentPtr);
    }
  }
}

std::vector<G4ReactionProduct> ExcitationHandler::ConvertResults(const G4FragmentVector& results) {
  std::vector<G4ReactionProduct> reactionProducts;
  reactionProducts.reserve(results.size());
  auto ionTable = G4ParticleTable::GetParticleTable()->GetIonTable();

  for (const auto& fragmentPtr : results) {
    auto fragmentDefinition = SpecialParticleDefinition(*fragmentPtr);
    if (fragmentDefinition == nullptr) {
      auto excitationEnergy = fragmentPtr->GetExcitationEnergy();
      auto level = fragmentPtr->GetFloatingLevelNumber();
      if (IsGroundState(*fragmentPtr)) {
        excitationEnergy = 0;
        level = 0;
      }
      fragmentDefinition = ionTable->GetIon(fragmentPtr->GetZ_asInt(), fragmentPtr->GetA_asInt(),
                                            excitationEnergy, G4Ions::FloatLevelBase(level));
    }
    // fragment wasn't found, ground state is created
    if (fragmentDefinition == nullptr) {
      fragmentDefinition = ionTable->GetIon(fragmentPtr->GetZ_asInt(), fragmentPtr->GetA_asInt(), 0, noFloat, 0);
      if (fragmentDefinition == nullptr) {
        throw std::runtime_error("ion table isn't created");
      }
      G4double ionMass = fragmentDefinition->GetPDGMass();
      if (fragmentPtr->GetMomentum().e() <= ionMass) {
        fragmentPtr->SetMomentum(G4LorentzVector(ionMass));
      } else {
        auto momentum = fragmentPtr->GetMomentum();
        G4double momentumModulus = std::sqrt((momentum.e() - ionMass) * (momentum.e() + ionMass));
        momentum.setVect(momentum.vect().unit() * momentumModulus);
        fragmentPtr->SetMomentum(momentum);
      }
    }

    reactionProducts.emplace_back(fragmentDefinition);
    reactionProducts.back().SetMomentum(fragmentPtr->GetMomentum().vect());
    reactionProducts.back().SetTotalEnergy((fragmentPtr->GetMomentum()).e());
    reactionProducts.back().SetFormationTime(fragmentPtr->GetCreationTime());
  }

  return reactionProducts;
}
