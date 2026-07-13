#ifndef DEEXCITATION_AAMCCHANDLERFACTORY_H_
#define DEEXCITATION_AAMCCHANDLERFACTORY_H_

#include <COLA.hh>

namespace cola {

  class AAMCCHandlerFactory final : public VConverterFactory {
   public:
    std::unique_ptr<VFilter> Create(const std::unordered_map<std::string, std::string>& meta_data) override;

    const std::string& GetFilterName() const override {
      static const std::string name{"AAMCCDeexcitationHandler"};
      return name;
    }
  };

}  // namespace cola

#endif  // DEEXCITATION_AAMCCHANDLERFACTORY_H_
