#include "Deexcitation/AAMCCHandlerFactory.h"

#include "Deexcitation/AAMCCHandlerConverter.h"
#include "Deexcitation/aamcc/include/DeexcitationHandler.hh"

#include <G4FermiDataTypes.hh>

#include <memory>
#include <optional>
#include <string>
#include <unordered_map>

namespace {

  double StodWithFactor(const std::string& value) {
    const double num = std::stod(value);
    if (value.ends_with("eV")) {
      return num * CLHEP::eV;
    }
    if (value.ends_with("MeV")) {
      return num * CLHEP::MeV;
    }
    if (value.ends_with("keV")) {
      return num * CLHEP::keV;
    }
    if (value.ends_with("GeV")) {
      return num * CLHEP::GeV;
    }
    return num;
  }

  struct Config {
    explicit Config(const std::unordered_map<std::string, std::string>& params) {
      if (const auto it = params.find("A"); it != params.end()) {
        atomic_mass = std::stoi(it->second);
      }
      if (const auto it = params.find("Z"); it != params.end()) {
        charge = std::stoi(it->second);
      }
      if (const auto it = params.find("lowerMfThreshold"); it != params.end()) {
        lower_mf_threshold = StodWithFactor(it->second);
      }
      if (const auto it = params.find("upperMfThreshold"); it != params.end()) {
        upper_mf_threshold = StodWithFactor(it->second);
      }
      if (const auto it = params.find("stableThreshold"); it != params.end()) {
        stable_threshold = StodWithFactor(it->second);
      }
      if (const auto it = params.find("minExForFermiBreakUp"); it != params.end()) {
        min_ex_for_fermi_break_up = StodWithFactor(it->second);
      }
    }

    std::optional<int> atomic_mass;
    std::optional<int> charge;
    std::optional<double> stable_threshold;
    std::optional<double> lower_mf_threshold;
    std::optional<double> upper_mf_threshold;
    std::optional<double> min_ex_for_fermi_break_up;
  };

}  // namespace

namespace cola {

  std::unique_ptr<VFilter> AAMCCHandlerFactory::Create(const std::unordered_map<std::string, std::string>& meta_data) {
    const Config config(meta_data);

    auto model = std::make_unique<DeexcitationHandler>();

    model->SetMaxAforFermiBreakUp(config.atomic_mass.value_or(static_cast<int>(MAX_A)));
    model->SetMaxZforFermiBreakUp(config.charge.value_or(static_cast<int>(MAX_Z)));
    model->SetExForMF(config.lower_mf_threshold.value_or(3 * CLHEP::MeV),
                      config.upper_mf_threshold.value_or(5 * CLHEP::MeV));

    if (config.stable_threshold.has_value()) {
      model->SetMinEx(*config.stable_threshold);
    }

    if (config.min_ex_for_fermi_break_up.has_value()) {
      model->SetMinExForFermiBreakUp(*config.min_ex_for_fermi_break_up);
    }

    return std::make_unique<AAMCCHandlerConverter>(std::move(model));
  }

}  // namespace cola
