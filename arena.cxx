/*
 * Sin-Yaw Wang <swang24@scu.edu>
 */
#include <bitset>
#include <vector>
#include <random>
#include <stdexcept>
#include "rng.h"
#include "dna.h"
#include "life.h"
#include "select.h"
#include "arena.h"

namespace csen79 {

/*
 * There are two way to keep scores.  We chose to do it individually.
 * Obviously, we could have picked random group of individuals and have them
 * "Compete" as a group.
 */
void Arena::compete() {
    for (auto& individual : population)
        for (int i = 0; i < EVENTPERROUND; i++) {
            competition.init();
            competition.play(individual); // Assuming play method modifies individual's score
        }
}

/*
 * Run after a round of competition, which can also be a "generation."
 * We select, arbitrarily, the upper half to preserve to the next generation.
 * Then we populate the rest with the offsprings of those winners.
 */
void Arena::evolve() {
    std::vector<Life> winners;
    
    // Select individuals with above-average scores
    long sumScore = 0;
    for (auto& individual : population)
        sumScore += individual.getScore();
    sumScore /= NPop;   // average score
    for (auto& individual : population)
        if (individual.getScore() >= sumScore)
            winners.push_back(individual);  // yes, make a copy

    population.clear();     // effectively wipe off the old generation
    // Add back winners, then their offsprings
    int i;
    for (i = 0; i < winners.size(); i++)
        population.push_back(winners[i]);
    for (; i < NPop; i++) {
        int p1 = rand_int(0, static_cast<int>(winners.size()) - 1);
        int p2 = rand_int(0, static_cast<int>(winners.size()) - 1);
        population.push_back(winners[p1].cross(winners[p2]));
    }
    // reset scores.  Let the next selection cycle begin.
    for (auto& individual : population)
        individual.newLife();
}

double Arena::matured() {
   double matches{0.0};
   for (auto& individual : population)
         matches += competition.matchRatio(individual);
    matches /= NPop;
    return matches;
}

}   // namespace csen79

