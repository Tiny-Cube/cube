#include "ukmc_io.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc, char** argv) {
  ukmc::require(argc == 4, "MODEL STEPS SEED required");
  std::ifstream input(argv[1]);
  auto model = ukmc::load(input);
  ukmc::Engine engine(*model, std::stoll(argv[3]));
  ukmc::Count checksum = 0;
  auto started = std::chrono::steady_clock::now();
  engine.run(std::stoll(argv[2]), ukmc::infinity,
             [&](const ukmc::Record& row) { checksum += row.event + 1; });
  double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
  std::cout << std::setprecision(17) << elapsed << ' ' << engine.steps << ' '
            << engine.time << ' ' << checksum;
  for (auto x : model->position()) std::cout << ' ' << x;
  for (auto x : model->state()) std::cout << ' ' << x;
  std::cout << '\n';
}
