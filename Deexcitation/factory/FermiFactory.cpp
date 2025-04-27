//
// Created by Artem Novikov on 12.06.2024.
//

#include <FermiBreakUp/Splitter.h>
#include <memory>
#include <optional>
#include <string>

#include <FermiBreakUp/FermiBreakUp.h>
#include <FermiBreakUp/nuclei_properties/data_storage/CSVNuclearMass.h>
#include <FermiBreakUp/nuclei_properties/NucleiProperties.h>

#include "Deexcitation/factory/FermiFactory.h"
#include "Deexcitation/converter/FermiConverter.h"
#include "FermiBreakUp/util/Cache.h"
#include "FermiBreakUp/util/DataTypes.h"

using namespace cola;

namespace {
  struct Config {
    std::optional<std::unique_ptr<fbu::FermiBreakUp::SplitCache>> cache;
    std::optional<std::string> nucleiCsv;
  };
}

cola::FermiConverter* FermiFactory::DoCreate(const std::map<std::string, std::string>& params) {
  Config config;
  if (auto it = params.find("cache"); it != params.end()) {
    const auto& [_, name] = *it;
    if (name == "simple") {
      config.cache = std::make_unique<fbu::SimpleCache<fbu::NucleiData, fbu::FragmentSplits>>();
    } else if (name == "lfu") {
      config.cache = std::make_unique<fbu::LFUCache<fbu::NucleiData, fbu::FragmentSplits>>(19 * 19 / 2);
    } else {
      throw std::runtime_error(R"(only "lfu" and "simple" caches are supported)");
    }
  }

  if (auto it = params.find("nucleiCsv"); it != params.end()) {
    const auto& [_, path] = *it;
    config.nucleiCsv = path;
  }

  auto model = config.cache.has_value() ? std::make_unique<fbu::FermiBreakUp>(*std::move(config.cache)) : std::make_unique<fbu::FermiBreakUp>();

  if (config.nucleiCsv.has_value()) {
    fbu::NucleiProperties::Reset(fbu::CSVNuclearMass(config.nucleiCsv.value()));
  }

  return new FermiConverter(std::move(model));
}
