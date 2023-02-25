#ifndef H_TIMER_STRUCT_H
#define H_TIMER_STRUCT_H

#include <chrono>

class Timer
{
private:
    std::chrono::time_point<std::chrono::steady_clock> start_time_;
    std::chrono::time_point<std::chrono::steady_clock> end_time_;
    std::chrono::duration<double, std::milli> elapsed_time_ = std::chrono::steady_clock::duration::zero();

public:
    void start() {
        start_time_ = std::chrono::steady_clock::now();
    }
    
    void stop() {
        end_time_ = std::chrono::steady_clock::now();
        elapsed_time_ += std::chrono::duration<double, std::milli>(end_time_ - start_time_);
    }

    double elapsed_time() const {
        return std::chrono::duration<double>(elapsed_time_).count();
    }
    
    Timer() = default;
    ~Timer() = default;
};

#endif
