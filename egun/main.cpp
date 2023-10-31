#include "olptm/core.hxx"

#include "constants.h"
#include <iostream>

using namespace lptm;

int main() {
  system_t system;

  body_t filaments[def::n_filamets];
  body_t aluminas[def::n_filamets];
  for (u32 i = 0; i < def::n_filamets; i++) {
    filaments[i] = body(def::fil_masses[i % 2], SpecificHeat::Tungsten,
                        Conductivity::Tungsten, 2000);
    aluminas[i] = body(def::alu_masses[i % 2], SpecificHeat::Alumina,
                       Conductivity::Alumina, 1000);
    add_body(system, filaments[i]);
    add_body(system, aluminas[i]);
    add_exchange(system, radiation(filaments[i], aluminas[i], 1.0,
                                   def::alumina_fialemnt[i % 2]));
  }

  f64 dt = 0.001;            // s
  const f64 sim_time = 10.0; // s
  f64 time = 0;              // s
  u32 iteration = 0;
  while (time < sim_time) {
    evaluate(system, dt);
    time += dt;
  }
  for (u32 i = 0; i < def::n_filamets; i++) {
    std::cout << system.temperatures[filaments[i].id] << ", ";
  }
  std::cout << std::endl;
  for (u32 i = 0; i < def::n_filamets; i++) {
    std::cout << system.temperatures[aluminas[i].id] << ", ";
  }
  std::cout << std::endl;
  return 0;
}
