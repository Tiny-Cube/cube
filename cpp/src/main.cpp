#include "ukmc_io.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc, char **argv) {
  try {
    ukmc::require(argc >= 2 && argc <= 5,
                  "usage: ukmc-cpp MODEL [STEPS [SEED [UNTIL]]]");
    auto integer = [](const char *s) {
      std::istringstream in(s);
      ukmc::Tokens t(in);
      auto v = t.integer();
      t.finish();
      return v;
    };
    auto number = [](const char *s) {
      std::istringstream in(s);
      ukmc::Tokens t(in);
      auto v = t.number();
      t.finish();
      return v;
    };
    const auto steps = argc > 2 ? integer(argv[2]) : 1000;
    const auto seed = argc > 3 ? integer(argv[3]) : 1;
    const auto until = argc > 4 ? number(argv[4]) : ukmc::infinity;
    ukmc::require(until >= 0 && (argc <= 4 || std::isfinite(until)),
                  "invalid time limit");
    std::ifstream input(argv[1]);
    ukmc::require(bool(input), "cannot open model");
    auto model = ukmc::load(input);
    ukmc::Engine engine(*model, seed);
    std::cout << "step,time,event,x,y,z";
    for (std::size_t i = 0; i < model->state().size(); ++i)
      std::cout << ",state" << i;
    std::cout << '\n' << std::setprecision(17);
    engine.run(steps, until, [](const ukmc::Record &r) {
      std::cout << r.step << ',' << r.time << ',' << r.event;
      for (auto x : r.position)
        std::cout << ',' << x;
      for (auto n : r.state)
        std::cout << ',' << n;
      std::cout << '\n';
    });
    std::cerr << "status=" << engine.status << " events=" << engine.steps
              << " time=" << std::setprecision(17) << engine.time << '\n';
  } catch (const std::exception &e) {
    std::cerr << "ukmc: " << e.what() << '\n';
    return 2;
  }
}
