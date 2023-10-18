#include "olptm/core.hxx"
#include <eigen3/Eigen/Core>
#include <fstream>
#include <iostream>

// formatter to save vector to CSV
const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                       Eigen::DontAlignCols, "\n", ", ");

using namespace lptm;

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cout << "Missing CSV file path\n";
    return 1;
  }
  body_t heat_src = create(180);
  body_t body_1 = create(0.05, 900, 23);
  body_t body_2 = create(0.05, 450, 23);
  body_t heat_snk = create(20);

  system_t system;

  add_body(system, heat_src);
  add_body(system, body_1);
  add_body(system, body_2);
  add_body(system, heat_snk);

  add_exchange(system, conduction(heat_src, body_1, Req(5.0 / 4.0, 400),
                                  Req(5.0 / 4.0, 200)));
  add_exchange(system, conduction(body_1, body_2, Req(5.0 / 4.0, 200),
                                  Req(5.0 / 4.0, 60)));
  add_exchange(system, conduction(body_2, heat_snk, Req(5.0 / 4.0, 60),
                                  Req(5.0 / 4.0, 400)));

  initialize(system);

  const f64 dt = 0.001;
  const u32 N = 10000;
  std::ofstream file(argv[1]);
  file << "#dt: " << dt << "\n";
  file << system.temperatures.format(CSVFormat) << "\n";
  for (u32 i = 0; i < N; i++) {
    evaluate(system, dt);
    file << system.temperatures.format(CSVFormat) << "\n";
  }
  std::cout << "Temperature heat source: " << system.temperatures[heat_src.id]
            << std::endl;
  std::cout << "Temperature heat sink: " << system.temperatures[heat_snk.id]
            << std::endl;
  std::cout << "Temperature b1: " << system.temperatures[body_1.id]
            << std::endl;
  std::cout << "Temperature b2: " << system.temperatures[body_2.id]
            << std::endl;
  return 0;
}
