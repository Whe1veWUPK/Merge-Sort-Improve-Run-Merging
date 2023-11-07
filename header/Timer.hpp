#pragma once
#include<chrono>
class Timer {
private:
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
public:
    void startTimer(); // start timer
    void endTimer(); // end timer
    void calculateTime(); //calculate the time
};