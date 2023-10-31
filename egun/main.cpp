#include "constants.h"

#include "olptm/core.hxx"
#include <chrono>
#include <fstream>
#include <iostream>

// formatter to save vector to CSV
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                       Eigen::DontAlignCols, "\n", ", ");

using namespace lptm;

exchange_t thermoionic_cooling(body_t emitter, f64 &Ik) {
  return [emitter, &Ik](system_t &sys) -> void {
      // formula to compute dE from Ik and emitter temeperature
  };
}
exchange_t filament_heating(body_t filament, f64 &Pf) {
  return
      [filament, &Pf](system_t &sys) -> void { sys.heats[filament.id] += Pf; };
}

// template <u32 n>
// exchange_t matrix_radiation(std::vector<body_t> bodies,
//                             const Eigen::Matrix<f64, n, n> equivR) {
//   u32 index[n];
//   for (u32 i = 0; i < n; i++) {
//     index[i] = bodies[i].id;
//   }
//   return [index, equivR](system_t &sys) -> void {
//     Eigen::Vector<f64, n> temp;
//     for (u32 i = 0; i < n; i++) {
//       temp[i] = _P4_(sys.temperatures[index[i]]);
//     }
//     Eigen::Vector<f64, n> heat = equivR * temp;
//     for (u32 i = 0; i < n; i++) {
//       sys.heats[index[i]] += heat(i);
//     }
//   };
// }

// in the arguments list should passed the optimized parameters and external values
bool initialise(system_t &system, f64 &Pf) {

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
    add_exchange(system, filament_heating(filaments[i], Pf));
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
  inner_cavity.push_back(emitter[0]);
  inner_cavity.push_back(upper_shields[0]);
  outer_cavity.push_back(emitter[0]);
  outer_cavity.push_back(emitter[0]);

  for (u32 i = 0; i < def::InnerCavity::n_bodies; i++) {
    for (u32 j = i + 1; j < def::InnerCavity::n_bodies; j++) {
      add_exchange(system,
                   radiation(inner_cavity[i], inner_cavity[j], 1,
                             def::InnerCavity::equivalent_resistance(i, j)));
    }
  }

  for (u32 i = 0; i < def::OuterCavity::n_bodies; i++) {
    for (u32 j = i + 1; j < def::OuterCavity::n_bodies; j++) {
      add_exchange(system,
                   radiation(outer_cavity[i], outer_cavity[j], 1,
                             def::OuterCavity::equivalent_resistance(i, j)));
    }
  }
  return true;
}


int main(int argc, char **argv) {
  system_t system;
  f64 Pf = 200;
  initialise(system, Pf);
  /*********************/
  /* SIMULATION CONFIG */
  /*********************/
  f64 dt = 0.001;            // s
  const f64 dT = 0.1;        // K
  const f64 sim_time = 10.0; // s
  bool csv_enabled = argc > 1;
  std::ofstream file;
  /*******************/
  /*   SYSTEM INFO   */
  /*******************/
  std::cout << "N Bodies: " << system.bodies.size() << "\n";
  std::cout << "N Exchanges: " << system.exchanges.size() << "\n";
  std::cout << "Desired dT: " << dT << "K\n";
  std::cout << "Simulation length: " << sim_time << "s\n";

  /*******************/
  /* SIMULATE SYSTEM */
  /*******************/
  f64 time = 0; // s
  u32 iteration = 0;
  if (csv_enabled) {
    file.open(argv[argc - 1]);
    file << iteration << ", " << time << ", "
         << system.temperatures.format(CSVFormat) << "\n";
  }
  std::chrono::time_point start = std::chrono::high_resolution_clock::now();
  while (time < sim_time) {
    evaluate(system, dT, dt, 0.0001, 1.0); // dT = 0.1 K
    time += dt;
    iteration++;
    if (csv_enabled) {
      file << iteration << ", " << time << ", "
           << system.temperatures.format(CSVFormat) << "\n";
    }
  }
  std::chrono::time_point stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Number of iterations: " << iteration << std::endl;
  std::cout << "Execution time: " << duration.count() << " us ("
            << duration.count() / iteration << " us x iter)\n";
  return 0;
}
