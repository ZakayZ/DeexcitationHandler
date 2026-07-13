#include "Deexcitation/AAMCCHandlerFactory.h"
#include "Deexcitation/G4HandlerFactory.h"

#include <CLHEP/Units/SystemOfUnits.h>
#include <COLA.hh>
#include <gtest/gtest.h>

#include <cmath>
#include <sstream>

namespace {

  constexpr auto kG4PipelineConfigXml = R"(<?xml version="1.0" encoding="UTF-8" ?>
<program>
    <generator name="generator"/>
    <converter name="G4DeexcitationHandler"/>
    <writer name="writer"/>
</program>
)";

  constexpr auto kAamccPipelineConfigXml = R"(<?xml version="1.0" encoding="UTF-8" ?>
<program>
    <generator name="generator"/>
    <converter name="AAMCCDeexcitationHandler"/>
    <writer name="writer"/>
</program>
)";

  class TestGenerator final : public cola::VGenerator {
   public:
    explicit TestGenerator(cola::EventParticles particles) : particles_(std::move(particles)) {}

    std::unique_ptr<cola::EventData> operator()() override {
      auto data = std::make_unique<cola::EventData>();
      data->particles = particles_;
      return data;
    }

   private:
    cola::EventParticles particles_;
  };

  class TestGeneratorFactory final : public cola::VGeneratorFactory {
   public:
    explicit TestGeneratorFactory(cola::EventParticles particles) : particles_(std::move(particles)) {}

    std::unique_ptr<cola::VFilter> Create(const std::unordered_map<std::string, std::string>& /*meta_data*/) override {
      return std::make_unique<TestGenerator>(particles_);
    }

    const std::string& GetFilterName() const override {
      static const std::string name{"generator"};
      return name;
    }

   private:
    cola::EventParticles particles_;
  };

  class TestWriter final : public cola::VWriter {
   public:
    explicit TestWriter(std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink) : sink_(std::move(sink)) {}

    void operator()(std::unique_ptr<cola::EventData>&& data) override { sink_->emplace_back(std::move(data)); }

   private:
    std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink_;
  };

  class TestWriterFactory final : public cola::VWriterFactory {
   public:
    explicit TestWriterFactory(std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink)
        : sink_(std::move(sink)) {}

    std::unique_ptr<cola::VFilter> Create(const std::unordered_map<std::string, std::string>& /*meta_data*/) override {
      return std::make_unique<TestWriter>(sink_);
    }

    const std::string& GetFilterName() const override {
      static const std::string name{"writer"};
      return name;
    }

   private:
    std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink_;
  };

  void RegisterTestPipeline(cola::MetaProcessor& meta_processor,
                            std::shared_ptr<std::vector<std::unique_ptr<cola::EventData>>> sink,
                            cola::EventParticles particles,
                            std::unique_ptr<cola::VConverterFactory> converter_factory) {
    meta_processor.Register(std::make_unique<TestGeneratorFactory>(std::move(particles)));
    meta_processor.Register(std::move(converter_factory));
    meta_processor.Register(std::make_unique<TestWriterFactory>(std::move(sink)));
  }

  void RunDeexcitationTest(const std::string& pipeline_xml, std::unique_ptr<cola::VConverterFactory> converter_factory,
                           const cola::Particle& particle) {
    auto sink = std::make_shared<std::vector<std::unique_ptr<cola::EventData>>>();
    cola::MetaProcessor meta_processor;
    RegisterTestPipeline(meta_processor, sink, {particle}, std::move(converter_factory));
    std::istringstream xml(pipeline_xml);
    cola::ColaRunManager manager(meta_processor.Parse(xml));
    manager.Run(1);

    ASSERT_EQ(sink->size(), 1u);
    EXPECT_GE((*sink)[0]->particles.size(), 1u);
  }

}  // namespace

TEST(TestModule, TestG4Handler) {
  const cola::Particle light_fragment{
      .position = cola::LorentzVector{},
      .momentum =
          cola::LorentzVector{
              .e = std::sqrt(3 * std::pow(100 * CLHEP::MeV, 2) + std::pow(4 * 938 * CLHEP::MeV, 2)),
              .x = 100 * CLHEP::MeV,
              .y = 100 * CLHEP::MeV,
              .z = 100 * CLHEP::MeV,
          },
      .pdg_code = cola::AZToPdg({4, 2}),
      .p_class = cola::ParticleClass::kSpectatorA,
  };

  const cola::Particle heavy_fragment{
      .position = cola::LorentzVector{},
      .momentum =
          cola::LorentzVector{
              .e = std::sqrt(3 * std::pow(100 * CLHEP::MeV, 2) + std::pow(5 * 938 * CLHEP::MeV, 2)),
              .x = 100 * CLHEP::MeV,
              .y = 100 * CLHEP::MeV,
              .z = 100 * CLHEP::MeV,
          },
      .pdg_code = cola::AZToPdg({5, 3}),
      .p_class = cola::ParticleClass::kSpectatorA,
  };

  RunDeexcitationTest(kG4PipelineConfigXml, std::make_unique<cola::G4HandlerFactory>(), light_fragment);
  RunDeexcitationTest(kG4PipelineConfigXml, std::make_unique<cola::G4HandlerFactory>(), heavy_fragment);
}

TEST(TestModule, TestAamccHandler) {
  const cola::Particle light_fragment{
      .position = cola::LorentzVector{},
      .momentum =
          cola::LorentzVector{
              .e = std::sqrt(3 * std::pow(100 * CLHEP::MeV, 2) + std::pow(4 * 938 * CLHEP::MeV, 2)),
              .x = 100 * CLHEP::MeV,
              .y = 100 * CLHEP::MeV,
              .z = 100 * CLHEP::MeV,
          },
      .pdg_code = cola::AZToPdg({4, 2}),
      .p_class = cola::ParticleClass::kSpectatorA,
  };

  const cola::Particle heavy_fragment{
      .position = cola::LorentzVector{},
      .momentum =
          cola::LorentzVector{
              .e = std::sqrt(3 * std::pow(100 * CLHEP::MeV, 2) + std::pow(5 * 938 * CLHEP::MeV, 2)),
              .x = 100 * CLHEP::MeV,
              .y = 100 * CLHEP::MeV,
              .z = 100 * CLHEP::MeV,
          },
      .pdg_code = cola::AZToPdg({5, 3}),
      .p_class = cola::ParticleClass::kSpectatorA,
  };

  RunDeexcitationTest(kAamccPipelineConfigXml, std::make_unique<cola::AAMCCHandlerFactory>(), light_fragment);
  RunDeexcitationTest(kAamccPipelineConfigXml, std::make_unique<cola::AAMCCHandlerFactory>(), heavy_fragment);
}

TEST(ModuleExport, LoadCOLAModuleExposesDeexcitationFactories) {
  auto module = std::unique_ptr<cola::VModule>(LoadCOLAModule());
  ASSERT_NE(module, nullptr);

  const auto filters = module->GetModuleFilters();
  ASSERT_EQ(filters.size(), 2u);
  ASSERT_TRUE(filters.contains("G4DeexcitationHandler"));
  ASSERT_TRUE(filters.contains("AAMCCDeexcitationHandler"));

  const auto* g4_factory = filters.at("G4DeexcitationHandler").get();
  ASSERT_NE(g4_factory, nullptr);
  EXPECT_EQ(g4_factory->GetFilterName(), "G4DeexcitationHandler");
  EXPECT_EQ(g4_factory->GetFilterType(), cola::FilterType::kConverter);

  const auto* aamcc_factory = filters.at("AAMCCDeexcitationHandler").get();
  ASSERT_NE(aamcc_factory, nullptr);
  EXPECT_EQ(aamcc_factory->GetFilterName(), "AAMCCDeexcitationHandler");
  EXPECT_EQ(aamcc_factory->GetFilterType(), cola::FilterType::kConverter);
}
