#include "ukmc_events.hpp"
#include <deque>
#include <iostream>
struct Queue {
  std::deque<std::string> jobs;
  ukmc::Count completed = 0, next_id = 0;
};
int main() {
  ukmc::EventModel<Queue> model(
      Queue{},
      [](const Queue &) {
        return std::vector<ukmc::Event<Queue>>{
            {[](const Queue &) { return 2.0; },
             [](Queue &q) {
               q.jobs.push_back("job-" + std::to_string(q.next_id++));
             },
             [](const Queue &) { return true; }, "arrival"},
            {[](const Queue &) { return 3.0; },
             [](Queue &q) {
               q.jobs.pop_front();
               ++q.completed;
             },
             [](const Queue &q) { return !q.jobs.empty(); }, "service"}};
      },
      [](const Queue &q) {
        return std::vector<ukmc::Count>{static_cast<ukmc::Count>(q.jobs.size()),
                                        q.completed};
      });
  ukmc::Engine engine(model, 42);
  engine.run(10000, 100, [](const ukmc::Record &) {});
  std::cout << model.data.jobs.size() << ' ' << model.data.completed << ' '
            << engine.steps << '\n';
}
