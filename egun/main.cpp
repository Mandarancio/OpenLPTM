#include "olptm/core.hxx"
#include <eigen3/Eigen/src/Core/Matrix.h>

using namespace lptm;

namespace Density {};

namespace SH {
const f64 Alumina = 0.1;
const f64 Tungsten = 0.1;
}; // namespace SH

namespace TC {
const f64 Alumina = 100;
const f64 Tungsten = 100;
}; // namespace TC

const u32 n_filamets = 8;

const Eigen::VectorX<f64> filaments_mass = {};
const Eigen::VectorX<f64> aluminas_mass = {};
const Eigen::VectorX<f64> alumina_fialemnt = {};

int main() {
  system_t system;

  body_t filaments[n_filamets];
  body_t aluminas[n_filamets];
  for (u32 i = 0; i < n_filamets; i++) {
    filaments[i] = body(filaments_mass[i], SH::Tungsten, TC::Tungsten, 2000);
    aluminas[i] = body(aluminas_mass[i], SH::Alumina, TC::Alumina, 1000);
    add_body(system, filaments[i]);
    add_body(system, aluminas[i]);
    add_exchange(
        system, radiation(filaments[i], aluminas[i], 0.5, alumina_fialemnt[i]));
  }

  return 0;
}
