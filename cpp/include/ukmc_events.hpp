#pragma once
#include "ukmc.hpp"
namespace ukmc {
// State is unconstrained: containers, user structs, graphs, continuous
// attributes, etc.
template <class State> struct Event {
  std::function<double(const State &)> rate;
  std::function<void(State &)> action;
  std::function<bool(const State &)> enabled = [](const State &) {
    return true;
  };
  std::string name;
};
template <class State> class EventModel : public Model {
  std::function<std::vector<Event<State>>(const State &)> factory_;
  std::function<std::vector<Count>(const State &)> encode_;
  mutable std::vector<Event<State>> events_;
  mutable std::vector<double> rates_;

public:
  State data;
  EventModel(
      State initial,
      std::function<std::vector<Event<State>>(const State &)> factory,
      std::function<std::vector<Count>(const State &)> encode =
          [](const State &) { return std::vector<Count>{}; })
      : factory_(std::move(factory)), encode_(std::move(encode)),
        data(std::move(initial)) {}
  std::vector<double> propensities() const override {
    auto events = factory_(data);
    std::vector<double> rates;
    for (const auto &e : events)
      rates.push_back(e.enabled(data) ? nonnegative(e.rate(data)) : 0);
    events_ = std::move(events);
    rates_ = rates;
    return rates;
  }
  void fire(std::size_t event) override {
    require(event < rates_.size() && rates_[event] > 0,
            "inactive or unevaluated event");
    events_[event].action(data);
    events_.clear();
    rates_.clear();
  }
  std::vector<Count> state() const override { return encode_(data); }
};
} // namespace ukmc
