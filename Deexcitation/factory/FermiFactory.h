#pragma once

#include <COLA.hh>

#include "Deexcitation/converter/FermiConverter.h"

namespace cola {
  class FermiFactory final : public cola::VFactory {
  public:
    cola::VFilter* create(const std::map<std::string, std::string>& params) final {
      return DoCreate(params);
    }

  private:
    cola::FermiConverter* DoCreate(const std::map<std::string, std::string>& params);
  };
} // namespace cola
