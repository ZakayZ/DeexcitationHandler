#include <FermiBreakUp/FermiBreakUp.h>

using namespace fbu;

int main() {
    auto model = FermiBreakUp();
    Particle p(10_m, 6_c, LorentzVector(1, 2, 3, 4));
    model.BreakItUp(p);
}
