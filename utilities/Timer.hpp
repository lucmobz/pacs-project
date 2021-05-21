#pragma once

#include <chrono>
#include <iostream>
#include <string>

namespace utl {
/// \brief Timer that gives the time elapsed when it exits a scope (use stop() for a normal timer).
struct Timer {
  using clock = std::chrono::high_resolution_clock;
  Timer() : _id(_global_id++), _start(clock::now()) {}
  Timer(const std::string& message) : _message(message), _start(clock::now()) {}
  ~Timer() {
    if (_do_display_at_destruction) {
      auto stop = clock::now();
      std::chrono::duration<double> time = stop - _start;
      std::clog << "time elapsed: " << time.count() << " ("
                << (_message.empty() ? std::to_string(_id) : _message) << ")\n";
    }
  }
  void stop() {
    auto stop = clock::now();
    std::chrono::duration<double> time = stop - _start;
    std::clog << "time elapsed: " << time.count() << " ("
              << (_message.empty() ? std::to_string(_id) : _message) << ")\n";
    _do_display_at_destruction = false;
  }

 private:
  inline static int _global_id{0};
  int _id{0};
  std::string _message;
  std::chrono::time_point<clock> _start;
  bool _do_display_at_destruction = true;
};
}