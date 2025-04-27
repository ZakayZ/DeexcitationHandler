#include <COLA.hh>
#include <CLHEP/Units/PhysicalConstants.h>

#include "Deexcitation/handler/ExcitationHandler.h"
#include "Deexcitation/DeexcitationModule.h"

class TestGenerator: public cola::VGenerator {
public:
  TestGenerator(const cola::EventParticles& particles) : particles_(particles) {}

  std::unique_ptr<cola::EventData> operator()() override {
    auto data = std::make_unique<cola::EventData>();
    data->particles = particles_;

    return data;
  }

private:
  const cola::EventParticles particles_;
};

class TestGeneratorFactory: public cola::VFactory {
public:
  cola::VFilter* create(const std::map<std::string, std::string>&) override {
    return new TestGenerator(particles);
  }

  cola::EventParticles particles;
};

class TestWriter: public cola::VWriter {
public:
  TestWriter(std::vector<std::unique_ptr<cola::EventData>>& holder) : events_(holder) {}

  void operator()(std::unique_ptr<cola::EventData>&& event) {
    events_.emplace_back(std::move(event));
  }

private:
  std::vector<std::unique_ptr<cola::EventData>>& events_;
};

class TestWriterFactory: public cola::VFactory {
public:
  cola::VFilter* create(const std::map<std::string, std::string>&) override {
    return new TestWriter(events);
  }

  std::vector<std::unique_ptr<cola::EventData>> events;
};

int main() {
  cola::MetaProcessor metaProcessor;
  TestGeneratorFactory* genFactory = new TestGeneratorFactory();
  TestWriterFactory* writerFactory = new TestWriterFactory();
  auto converter = new cola::FermiFactory();
  metaProcessor.reg(std::unique_ptr<cola::VFactory>(genFactory), "generator", cola::FilterType::generator);
  metaProcessor.reg(std::unique_ptr<cola::VFactory>(converter), "converter", cola::FilterType::converter);
  metaProcessor.reg(std::unique_ptr<cola::VFactory>(writerFactory), "writer", cola::FilterType::writer);
  auto manager = cola::ColaRunManager(metaProcessor.parse("config.xml"));
  manager.run(1);
}
