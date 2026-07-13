#ifndef DEEXCITATION_AAMCCHANDLERCONVERTER_H_
#define DEEXCITATION_AAMCCHANDLERCONVERTER_H_

#include <COLA.hh>

#include <memory>

class DeexcitationHandler;

namespace cola {
  class AAMCCHandlerConverter final : public cola::VConverter {
   public:
    explicit AAMCCHandlerConverter(std::unique_ptr<DeexcitationHandler>&& model);

    std::unique_ptr<cola::EventData> operator()(std::unique_ptr<cola::EventData>&& data) final;

   private:
    std::unique_ptr<DeexcitationHandler> model_;
  };
}  // namespace cola

#endif  // DEEXCITATION_AAMCCHANDLERCONVERTER_H_
