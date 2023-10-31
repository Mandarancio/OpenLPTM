#include "olptm/core.hxx"

#include "constants.h"
#include <iostream>

using namespace lptm;

exchange_t thermoionic_cooling(body_t &emitter, f64 &Ik) {
  return [&emitter, &Ik](system_t &sys) -> void {
    // formula to compute dE from Ik and emitter temeperature
  };
}
exchange_t filament_heating(body_t &filament, f64 &Pf) {
  return
      [&filament, &Pf](system_t &sys) -> void { sys.heats[filament.id] += Pf; };
}

int main() {
  const f64 dt = 0.001;      // s
  const f64 sim_time = 10.0; // s

  system_t system;

  body_t filaments[def::n_filamets];
  body_t aluminas[def::n_filamets];
  body_t lower_shields[def::n_lower_shields];
  body_t upper_shields[def::n_upper_shields];
  body_t oil, nosecone, prolongator;

  std::vector<body_t> inner_cavity;
  std::vector<body_t> outer_cavity;

  // initialisation filaments and aluminas bodies and exchanges
  for (u32 i = 0; i < def::n_filamets; i++) {
    filaments[i] = body(def::fil_masses[i % 2], SpecificHeat::Tungsten,
                        Conductivity::Tungsten, 2000);
    aluminas[i] = body(def::alu_masses[i % 2], SpecificHeat::Alumina,
                       Conductivity::Alumina, 1000);
    add_body(system, filaments[i]);
    add_body(system, aluminas[i]);
    add_exchange(system, radiation(filaments[i], aluminas[i], 1.0,
                                   def::alumina_fialemnt[i % 2]));
    inner_cavity.push_back(aluminas[i]);
  }

  // initialisation lower shields bodies
  for (u32 i = 0; i < def::n_lower_shields; i++) {
    lower_shields[i] = body(def::ls_masses[i], SpecificHeat::Molybdenum,
                            Conductivity::Molybdenum, 1000);
    add_body(system, lower_shields[i]);
    outer_cavity.push_back(lower_shields[i]);
  }
  
  // initialisation nosecone and prolongator
  nosecone = body(def::ns_mass, SpecificHeat::Molybdenum,
                  Conductivity::Molybdenum, 1000);
  prolongator = body(def::pr_mass, SpecificHeat::Molybdenum,
                     Conductivity::Molybdenum, 1000);
  add_body(system, nosecone);
  add_body(system, prolongator);
  outer_cavity.push_back(nosecone);
  outer_cavity.push_back(prolongator);
  
  // initialisation upper shields bodies
  for (u32 i = 0; i < def::n_upper_shields; i++) {
    upper_shields[i] = body(def::us_masses[i], SpecificHeat::Molybdenum,
                            Conductivity::Molybdenum, 1000);
    add_body(system, upper_shields[i]);
    outer_cavity.push_back(upper_shields[i]);
  }

  // initialise oil
  oil = constant_temperature_body(900, Conductivity::Oil);
  add_body(system, oil);

  // add oil - nosecone / prolongator exchanges
  add_exchange(system, radiation(nosecone, oil, 1, def::oil_nosecone));
  add_exchange(system, radiation(prolongator, oil, 1, def::oil_prolongator));

  // initialisation emitter bodies
  body_t emitter[def::n_em_segments];
  for (u32 i = 0; i < def::n_em_segments; i++) {
    emitter[i] = body(def::em_masses[i], SpecificHeat::ImprTungsten,
                      Conductivity::ImprTungsten, 1200);
    add_body(system, emitter[i]);
  }
  outer_cavity.push_back(emitter[0]);
  outer_cavity.push_back(emitter[0]);

  /*******************/
  /*   SYSTEM INFO   */
  /*******************/
  std::cout << "N Bodies: " << system.bodies.size() << "\n";
  std::cout << "N Exchanges: " << system.exchanges.size() << "\n";
  std::cout << "Outer cavity N Bodies: "<< outer_cavity.size() << "\n";
  std::cout << "Inner cavity N Bodies: "<< inner_cavity.size() << "\n";
  std::cout << "Static dt: " << dt << "s\n";
  std::cout << "Simulation length: " << sim_time << "s\n";

  /*******************/
  /* SIMULATE SYSTEM */
  /*******************/

  f64 time = 0; // s
  u32 iteration = 0;
  while (time < sim_time) {
    evaluate(system, dt);
    time += dt;
  }
  return 0;
}
