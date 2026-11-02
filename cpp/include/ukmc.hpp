#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ukmc {
using Count = std::int64_t;
constexpr Count max_count = 9007199254740991LL;
constexpr double infinity = std::numeric_limits<double>::infinity();
inline void require(bool ok, const char *message) {
  if (!ok)
    throw std::invalid_argument(message);
}
inline Count count(Count n) {
  require(n >= 0 && n <= max_count, "count outside [0, 2^53-1]");
  return n;
}
inline double nonnegative(double a) {
  require(std::isfinite(a) && a >= 0, "rate must be finite and nonnegative");
  return a;
}
struct PortableRNG {
  Count state;
  explicit PortableRNG(Count seed = 1)
      : state(count(seed) % 2147483646LL + 1) {}
  double uniform() {
    state = (48271LL * state) % 2147483647LL;
    return double(state) / 2147483647.0;
  }
};
struct Model {
  virtual ~Model() = default;
  virtual std::vector<double> propensities() const = 0;
  virtual void fire(std::size_t event) = 0;
  virtual std::vector<Count> state() const { return {}; }
  virtual std::array<double, 3> position() const { return {}; }
};
struct Record {
  Count step;
  double time;
  Count event;
  std::array<double, 3> position;
  std::vector<Count> state;
};
class Engine {
  Model &model_;
  PortableRNG rng_;
  std::function<double()> random_;
  double uniform() {
    const auto u = random_ ? random_() : rng_.uniform();
    require(std::isfinite(u) && u > 0 && u < 1,
            "RNG must return a value in (0,1)");
    return u;
  }

public:
  double time = 0;
  Count steps = 0;
  std::string status = "ready";
  explicit Engine(Model &model, Count seed = 1,
                  std::function<double()> random = {})
      : model_(model), rng_(seed), random_(std::move(random)) {}
  Record snapshot(Count event = -1) const {
    return {steps, time, event, model_.position(), model_.state()};
  }
  bool step(double until = infinity, Count *selected = nullptr) {
    require(!std::isnan(until) && until >= time,
            "until must be >= current time");
    if (time == until) {
      status = "time_limit";
      return false;
    }
    const auto rates = model_.propensities();
    double total = 0;
    for (auto a : rates)
      total += nonnegative(a);
    require(std::isfinite(total), "total propensity overflow");
    if (total == 0) {
      status = "absorbed";
      if (std::isfinite(until))
        time = until;
      return false;
    }
    const auto dt = -std::log(uniform()) / total;
    if (dt > until - time) {
      time = until;
      status = "time_limit";
      return false;
    }
    require(std::isfinite(dt) && std::isfinite(time + dt) && time + dt > time,
            "waiting time outside floating-point range");
    const auto threshold = uniform() * total;
    double cumulative = 0;
    std::size_t chosen = 0;
    for (std::size_t i = 0; i < rates.size(); ++i) {
      if (rates[i] > 0)
        chosen = i;
      cumulative += rates[i];
      if (threshold < cumulative)
        break;
    }
    model_.fire(chosen);
    time += dt;
    ++steps;
    status = "running";
    if (selected)
      *selected = static_cast<Count>(chosen);
    return true;
  }
  void run(Count max_steps, double until,
           const std::function<void(const Record &)> &observe) {
    count(max_steps);
    require(!std::isnan(until) && until >= time,
            "until must be >= current time");
    observe(snapshot());
    for (Count i = 0; i < max_steps; ++i) {
      const auto before = time;
      Count event = -1;
      if (!step(until, &event)) {
        if (time != before)
          observe(snapshot());
        return;
      }
      observe(snapshot(event));
    }
    status = time == until ? "time_limit" : "step_limit";
  }
};
struct Edge {
  Count source, target;
  double rate;
  std::array<double, 3> displacement{};
};
class Graph : public Model {
  std::vector<Edge> edges_;
  Count site_;
  std::array<double, 3> xyz_{};

public:
  Graph(Count sites, std::vector<Edge> edges, Count initial = 0)
      : edges_(std::move(edges)), site_(initial) {
    require(count(sites) > 0 && count(initial) < sites,
            "invalid graph size or initial site");
    for (const auto &e : edges_) {
      require(count(e.source) < sites && count(e.target) < sites,
              "edge endpoint out of bounds");
      nonnegative(e.rate);
      for (double x : e.displacement)
        require(std::isfinite(x), "non-finite displacement");
    }
  }
  std::vector<double> propensities() const override {
    std::vector<double> result;
    result.reserve(edges_.size());
    for (const auto &e : edges_)
      result.push_back(e.source == site_ ? e.rate : 0);
    return result;
  }
  void fire(std::size_t event) override {
    require(event < edges_.size(), "event out of bounds");
    const auto &e = edges_[event];
    require(e.source == site_ && e.rate > 0, "inactive event");
    auto next = xyz_;
    for (std::size_t k = 0; k < 3; ++k) {
      next[k] += e.displacement[k];
      require(std::isfinite(next[k]), "position overflow");
    }
    xyz_ = next;
    site_ = e.target;
  }
  std::vector<Count> state() const override { return {site_}; }
  std::array<double, 3> position() const override { return xyz_; }
};
struct Reaction {
  double rate;
  std::vector<Count> reactants, products;
};
class ReactionNetwork : public Model {
  std::vector<Count> populations_;
  std::vector<Reaction> reactions_;

public:
  ReactionNetwork(std::vector<Count> populations,
                  std::vector<Reaction> reactions)
      : populations_(std::move(populations)), reactions_(std::move(reactions)) {
    require(!populations_.empty(), "at least one species required");
    for (auto n : populations_)
      count(n);
    for (const auto &r : reactions_) {
      nonnegative(r.rate);
      require(r.reactants.size() == populations_.size() &&
                  r.products.size() == populations_.size(),
              "stoichiometry dimension mismatch");
      for (auto n : r.reactants)
        count(n);
      for (auto n : r.products)
        count(n);
    }
  }
  std::vector<double> propensities() const override {
    std::vector<double> result;
    for (const auto &r : reactions_) {
      bool enabled = r.rate > 0;
      for (std::size_t i = 0; i < populations_.size(); ++i)
        if (populations_[i] < r.reactants[i])
          enabled = false;
      double a = enabled ? r.rate : 0;
      if (enabled)
        for (std::size_t i = 0; i < populations_.size(); ++i) {
          const auto n = populations_[i],
                     order = std::min(r.reactants[i], n - r.reactants[i]);
          for (Count k = 1; k <= order; ++k)
            a = nonnegative(a * (double(n - k + 1) / double(k)));
        }
      result.push_back(a);
    }
    return result;
  }
  void fire(std::size_t event) override {
    require(event < reactions_.size(), "event out of bounds");
    const auto &r = reactions_[event];
    require(r.rate > 0, "inactive reaction");
    auto next = populations_;
    for (std::size_t i = 0; i < next.size(); ++i) {
      require(next[i] >= r.reactants[i], "inactive reaction");
      next[i] = count(next[i] - r.reactants[i] + r.products[i]);
    }
    populations_ = std::move(next);
  }
  std::vector<Count> state() const override { return populations_; }
  std::array<double, 3> position() const override { return {}; }
};
constexpr double kb = 8.617333262145e-5, hbar = 6.582119569e-16,
                 pi = 3.14159265358979323846;
inline double thermal(double t) {
  require(std::isfinite(t) && t > 0, "positive finite temperature required");
  return kb * t;
}
inline double marcus(double coupling, double lambda, double dg,
                     double t = 300) {
  const auto kt = thermal(t);
  require(std::isfinite(coupling) && std::isfinite(lambda) && lambda > 0 &&
              std::isfinite(dg),
          "invalid Marcus energies");
  if (coupling == 0)
    return 0;
  const auto x = (dg + lambda) / std::sqrt(lambda * kt);
  return nonnegative(
      std::exp(std::log(2 * pi / hbar) + 2 * std::log(std::abs(coupling)) -
               0.5 * std::log(4 * pi * lambda * kt) - x * x / 4));
}
inline double arrhenius(double attempt, double barrier, double t = 300) {
  return nonnegative(attempt) * std::exp(-nonnegative(barrier) / thermal(t));
}
inline double miller_abrahams(double attempt, double distance,
                              double localization, double dg, double t = 300) {
  const auto kt = thermal(t);
  nonnegative(attempt);
  nonnegative(distance);
  require(std::isfinite(localization) && localization > 0 && std::isfinite(dg),
          "invalid hopping parameters");
  return attempt *
         std::exp(-2 * distance / localization - std::max(dg, 0.0) / kt);
}
} // namespace ukmc
