#include "progress.hpp"
#include "iostream"
#include "string"
#include "math.h"
#include "thread"
#include "chrono"

progress::progBar::progBar(double min, double max) : min(min), max(max) {
    state = min;
    startTime = std::chrono::high_resolution_clock::now();
    std::string toPrint;
    std::string bar(barLen, '_');
    toPrint = "|" + bar + "| " + std::to_string(100*(state-min)/(max-min)) + "%";
    prevStr = toPrint;
    std::cout << "Progress:" << std::endl;
    std::cout << toPrint << std::flush;
}

void progress::progBar::update(double status) {
    state = status;
    clock now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
    double frac = (state-min)/(max-min);
    auto ETA = (double)duration/frac - (double)duration;
    int filled = std::floor(barLen*frac);
    std::string toPrint;
    if (filled > 0 && barLen-filled > 0) {
        std::string bar1(filled, '#');
        std::string bar2(barLen - filled, '_');

        toPrint = "|" + bar1 + bar2 + "| " + std::to_string(100 * frac) + "%; Elapsed Time:" + std::to_string(1e-3*(double)duration/60) + "min; ETA: " + std::to_string(1e-3*ETA/60) + "min";
        std::string clear((int) prevStr.length(), '\b');
        std::cout << clear;
        std::cout << toPrint << std::flush;
        prevStr = toPrint;
    }
}

void progress::progBar::add(double extra) {
    state += extra;
}

progress::progBar::progBar(int num) {
    min = 0;
    state = min;
    max = num;
//    mainThread = std::thread;
}

void progress::progBar::mainLoop() {
    std::string toPrint;
    std::string bar(barLen, '_');
    toPrint = "|" + bar + "| " + std::to_string(100 * (state - min) / (max - min)) + "%";
    prevStr = toPrint;
    std::cout << "Progress:" << std::endl;
    std::cout << toPrint << std::flush;

    cont = true;
    while (cont) {
        clock now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime).count();
        double frac = (state-min)/(max-min);
        auto ETA = (double)duration/frac - (double)duration;

        int filled = std::floor(barLen*frac);
        std::string toPrint;
        std::string bar1(filled, '#');
        std::string bar2(barLen - filled, '_');

        toPrint = "|" + bar1 + bar2 + "| " + std::to_string(100 * frac) + "%; Elapsed Time:" + std::to_string(1e-3*(double)duration/60) + "min; ETA: " + std::to_string(1e-3*ETA/60) + "min";
        std::string clear((int) prevStr.length(), '\b');
        std::cout << clear;
        std::cout << toPrint << std::flush;
        prevStr = toPrint;
        std::this_thread::sleep_for(std::chrono::milliseconds (100));
    }

}

std::thread progress::progBar::start() {
    startTime = std::chrono::high_resolution_clock::now();
    return std::thread(&progBar::mainLoop, this);
}

void progress::progBar::finish() {
    cont = false;
}

void progress::progBar::reset(int num) {
    min = 0;
    state = min;
    max = num;
}
