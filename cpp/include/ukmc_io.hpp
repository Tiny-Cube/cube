#pragma once
#include "ukmc.hpp"
#include <cstdlib>
#include <istream>
#include <regex>
#include <sstream>
namespace ukmc {
class Tokens {
  std::istringstream stream_;

public:
  explicit Tokens(std::istream &input) {
    std::string text, line;
    while (std::getline(input, line))
      text += line.substr(0, line.find('#')) + " ";
    stream_.str(text);
  }
  std::string next() {
    std::string token;
    require(bool(stream_ >> token), "truncated KMC1 input");
    return token;
  }
  Count integer() {
    const auto s = next();
    require(std::regex_match(s, std::regex("[0-9]+")),
            "expected unsigned decimal integer");
    return count(std::stoll(s));
  }
  double number() {
    const auto s = next();
    require(
        std::regex_match(
            s, std::regex(
                   R"([+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?)")),
        "expected decimal number");
    // strtod permits underflow to zero; validation rejects overflow downstream.
    return std::strtod(s.c_str(), nullptr);
  }
  void finish() {
    std::string extra;
    require(!(stream_ >> extra), "unexpected trailing tokens");
  }
};
inline std::unique_ptr<Model> load(std::istream &input) {
  Tokens t(input);
  require(t.next() == "KMC1", "expected KMC1 header");
  const auto kind = t.next();
  const auto n = t.integer(), m = t.integer();
  std::unique_ptr<Model> model;
  if (kind == "graph") {
    const auto initial = t.integer();
    std::vector<Edge> edges;
    for (Count j = 0; j < m; ++j) {
      const auto source = t.integer(), target = t.integer();
      const auto rate = t.number();
      const auto x = t.number(), y = t.number(), z = t.number();
      edges.push_back({source, target, rate, {x, y, z}});
    }
    model = std::make_unique<Graph>(n, std::move(edges), initial);
  } else if (kind == "reaction") {
    std::vector<Count> populations;
    for (Count i = 0; i < n; ++i)
      populations.push_back(t.integer());
    std::vector<Reaction> reactions;
    for (Count j = 0; j < m; ++j) {
      Reaction r{t.number(), {}, {}};
      for (Count i = 0; i < n; ++i)
        r.reactants.push_back(t.integer());
      for (Count i = 0; i < n; ++i)
        r.products.push_back(t.integer());
      reactions.push_back(std::move(r));
    }
    model = std::make_unique<ReactionNetwork>(std::move(populations),
                                              std::move(reactions));
  } else
    throw std::invalid_argument("unknown model kind");
  t.finish();
  return model;
}
} // namespace ukmc
