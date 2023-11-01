#ifndef QP_TIMER_H_
#define QP_TIMER_H_

#include <sys/time.h>

class timer {
 public:
  timer() {
    gettimeofday(&time,nullptr);
    time_start = 0.0;
    time_get = 0.0;
  }
  ~timer() {}
  void start() {
    gettimeofday(&time,nullptr);
    time_start = static_cast<double>(time.tv_sec)*1000.0 + 0.001*static_cast<double>(time.tv_usec);
  }
  double get() {
    gettimeofday(&time,nullptr);
    gettimeofday(&time,nullptr);
    time_get = static_cast<double>(time.tv_sec)*1000.0 + 0.001*static_cast<double>(time.tv_usec);
    return time_get - time_start;
  }
 private:
  double time_start;
  double time_get;
  timeval time;
};

#endif // QP_TIMER_H_