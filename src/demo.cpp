#include "olptm/constants.h"
#include "olptm/core.hxx"

#include <eigen3/Eigen/Core>

#include <chrono>
#include <fstream>
#include <iostream>

#define INV_CC(K) [](f64 T) -> f64 { return 1.0 / K; }
#define CTC(K) [](f64 T) -> f64 { return  K; }

// formatter to save vector to CSV
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                       Eigen::DontAlignCols, "\n", ", ");

using namespace lptm;

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Missing CSV file path\n";
    return 1;
  }

  bool static_dt = true;
  if (argc == 3) {
    if (strcmp(argv[1], "--dynamic") == 0)
      static_dt = false;
    else if (strcmp(argv[1], "--static") == 0)
      static_dt = true;
    else {
      std::cout << "mode can be only `--dynamic` or `--static`\n";
      return 1;
    }
  }

  body_t heat_src = constant_temperature_body(180, CC(ThermalConductivity::Cu));
  body_t body_1 =
      body(0.05, CC(SpecificHeat::Al), CC(ThermalConductivity::Al), 23);
  body_t body_2 =
      body(0.05, CC(SpecificHeat::Cu), CC(ThermalConductivity::Cu), 23);
  body_t heat_snk = constant_temperature_body(20, CC(ThermalConductivity::Cu));

  system_t system;

  add_body(system, heat_src);
  add_body(system, body_1);
  add_body(system, body_2);
  add_body(system, heat_snk);

  add_exchange(system, conduction(heat_src, body_1,
                                  Req(5.0 / 4.0, heat_src.thermal_conductivy),
                                  Req(5.0 / 4.0, body_1.thermal_conductivy)));
  add_exchange(system, conduction(body_1, body_2,
                                  Req(5.0 / 4.0, body_1.thermal_conductivy),
                                  Req(5.0 / 4.0, body_2.thermal_conductivy)));
  add_exchange(system, conduction(body_2, heat_snk,
                                  Req(5.0 / 4.0, body_1.thermal_conductivy),
                                  Req(5.0 / 4.0, heat_snk.thermal_conductivy)));

  f64 dt = 0.001;            // s
  const f64 sim_time = 10.0; // s
  f64 time = 0;              // s
  u32 iteration = 0;
  std::ofstream file(argv[argc - 1]);
  file << iteration << ", " << time << ", "
       << system.temperatures.format(CSVFormat) << "\n";
  std::chrono::time_point start = std::chrono::high_resolution_clock::now();
  while (time < sim_time) {
    if (static_dt) {
      evaluate(system, dt);
    } else {
      evaluate(system, 0.01, dt, 0.0001, 1.0);
    }
    time += dt;
    file << ++iteration << ", " << time << ", "
         << system.temperatures.format(CSVFormat) << "\n";
  }
  std::chrono::time_point stop = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Final time: " << time << " s\n";
  std::cout << "Number of iterations: " << iteration << std::endl;
  std::cout << "Execution time: " << duration.count() << " us ("
            << duration.count() / iteration << " us x iter)\n";
  std::cout << "Temperature b1: " << system.temperatures[body_1.id] << "C\n";
  std::cout << "Temperature b2: " << system.temperatures[body_2.id] << "C\n";
  return 0;
}
