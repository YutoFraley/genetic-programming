/*
 * Sin-Yaw Wang <swang24@scu.edu>
 */
#include <iostream>
#include <iomanip>
#include <random>
#include <stdexcept>
#include "rng.h"
#include "dna.h"
#include "life.h"
#include "select.h"
#include "arena.h"

// Global RNG implementation
// anonymous namespace is visible only to this file
namespace {
    std::mt19937& global_rng() {
        static std::mt19937 rng(std::random_device{}());
        return rng;
    }
}

using namespace csen79;

int main() {
    // RNG seeded automatically in rng.h via random_device
    Arena tribe;    // constructor populates them with random "genes."
    double avgScore;
    double maturity;
    do {
        maturity = tribe.matured();
        tribe.compete();
        avgScore = tribe.averageScore();
        int n = tribe.highGeneration();
        std::cout << "Generation: " << std::setw(3) << n << 
            ", Average score: " << std::setw(8) << std::fixed << std::setprecision(4) << avgScore << 
            ", Maturity: " << maturity << std::endl;
        tribe.evolve();
    } while (maturity <= Arena::MATURITY);
    // There needs a way to display the end population.
    return 0;
}

// RNG implementation
int rand_int(int lo, int hi) {
    if (hi <= lo) return lo;
    std::uniform_int_distribution<int> d(lo, hi);
    return d(global_rng());
}

bool rand_bool() {
    std::bernoulli_distribution d(0.5);
    return d(global_rng());
}
