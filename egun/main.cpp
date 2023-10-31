#include "olptm/core.hxx"

#include <eigen3/Eigen/Core>
#include <iostream>

using namespace lptm;

namespace Density {
const f64 Tungsten = 1.93e4;
const f64 Alumina = 1.93e4;
}; // namespace Density

namespace SH {
const f64 Alumina = 134.4;
const f64 Tungsten = 134.4;
}; // namespace SH

namespace TC {
const f64 Alumina = 170;
const f64 Tungsten = 170;
}; // namespace TC

namespace ME {
const f64 Tungsten = 0.5;
const f64 Alumina = 0.5;
}; // namespace ME

const u32 n_filamets = 8;

f64 radiationResistance(f64 e_b1, f64 A_b1, f64 e_b2, f64 A_b2) {
  double energy_resistance =
      (1 - e_b1) / (e_b1 * A_b1) + 1 / A_b1 + (1 - e_b2) / (e_b2 * A_b2);
  return 1. / energy_resistance;
}

const f64 fil_masses[2] = {
    Density::Tungsten * 1.03e-7,
    Density::Tungsten * 1.10e-7,
};
const f64 alu_masses[2] = {
    Density::Alumina * 8.190e-7,
    Density::Alumina * 8.758e-7,
};
const f64 alumina_fialemnt[2] = {
    radiationResistance(ME::Tungsten, 592.36e-6, ME::Alumina, 961.3e-6),
    radiationResistance(ME::Tungsten, 644.40e-6, ME::Alumina, 1027.9e-6),
};

int main() {
  system_t system;

  body_t filaments[n_filamets];
  body_t aluminas[n_filamets];
  for (u32 i = 0; i < n_filamets; i++) {
    filaments[i] = body(fil_masses[i % 2], SH::Tungsten, TC::Tungsten, 2000);
    aluminas[i] = body(alu_masses[i % 2], SH::Alumina, TC::Alumina, 1000);
    add_body(system, filaments[i]);
    add_body(system, aluminas[i]);
    add_exchange(system, radiation(filaments[i], aluminas[i], 1.0,
                                   alumina_fialemnt[i % 2]));
  }

  f64 dt = 0.001;            // s
  const f64 sim_time = 10.0; // s
  f64 time = 0;              // s
  u32 iteration = 0;
  while (time < sim_time) {
    evaluate(system, dt);
    time += dt;
  }
  for (u32 i = 0; i < n_filamets; i++) {
    std::cout << system.temperatures[filaments[i].id] << ", ";
  }
  std::cout << std::endl;
  for (u32 i = 0; i < n_filamets; i++) {
    std::cout << system.temperatures[aluminas[i].id] << ", ";
  }
  std::cout << std::endl;
  return 0;
}
