#include "ukmc_events.hpp"
#include "ukmc_io.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
void check(bool ok) {
  if (!ok)
    throw std::runtime_error("test failed");
}
int main() {
  using namespace ukmc;
  PortableRNG rng(0);
  check(rng.uniform() == 48271.0 / 2147483647.0);
  Graph g(1, {{0, 0, 2, {}}, {0, 0, 3, {}}});
  std::vector<double> u{std::exp(-1.0), .8, std::exp(-1.0)};
  std::size_t pos = 0;
  Engine engine(g, 1, [&] { return u.at(pos++); });
  std::vector<Record> rows;
  engine.run(10, .3, [&](const auto &r) { rows.push_back(r); });
  check(rows.size() == 3 && rows[1].event == 1 && rows[2].event == -1 &&
        engine.steps == 1 && engine.time == .3);
  ReactionNetwork reaction({10, 0},
                           {{.1, {2, 0}, {0, 1}}, {.3, {0, 1}, {2, 0}}});
  check(std::abs(reaction.propensities()[0] - 4.5) < 1e-12);
  Engine re(reaction, 99);
  re.run(2000, infinity,
         [](const auto &r) { check(r.state[0] + 2 * r.state[1] == 10); });
  bool rejected = false;
  try {
    Graph bad(1, {{0, 0, -1, {}}});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected);
  struct State {
    std::vector<std::string> tasks;
    double value;
  };
  EventModel<State> custom({{"a", "b", "c"}, 1.25}, [](const State &s) {
    std::vector<Event<State>> events;
    for (std::size_t i = 0; i < s.tasks.size(); ++i)
      events.push_back({[](const State &) { return 1.0; },
                        [i](State &d) {
                          d.tasks.erase(d.tasks.begin() + i);
                          d.value *= .5;
                        },
                        [](const State &) { return true; }, "complete"});
    return events;
  });
  Engine ce(custom, 13);
  ce.run(10, infinity, [](const auto &) {});
  check(ce.steps == 3 && ce.status == "absorbed" &&
        custom.data.value == .15625);
  struct Bare : Model {
    int n = 0;
    std::vector<double> propensities() const override {
      return n < 2 ? std::vector<double>{1} : std::vector<double>{};
    }
    void fire(std::size_t) override { ++n; }
  } bare;
  Engine be(bare);
  be.run(10, infinity, [](const auto &) {});
  check(bare.n == 2);
  check(std::abs(marcus(.01, .2, .05) / marcus(.01, .2, -.05) -
                 std::exp(-.05 / (kb * 300))) < 1e-13);
  std::cout << std::setprecision(17) << marcus(.01, .2, .05) << ' '
            << arrhenius(1e12, .1) << ' ' << miller_abrahams(1e12, 1, 2, .1)
            << '\n';
}
