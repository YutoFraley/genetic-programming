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

Select::Select(): favors{0x0f, 0xf0} {
    for (int c = 0; c < Life::NChromo; c++) 
        maxMatches[c] = favors[c].count();
}

/*
 * new score for a round of play for an individual
 * First see how much the gene matches the preferred
 * Compute the score favoring the matched gene, add some luck element
 */
void Select::play(Life& player) {
    double ratio{matchRatio(player)};
    long luck = rand_int(-100, 100);
    // Simulate playing the game and update player's score
    player.newScore((ratio - 0.5) * 30 + luck);
}

// % of gene codes that matched our preference
double Select::matchRatio(const Life& player) {
    double ratio{0.0};
    for (int c = 0; c < Life::NChromo; c++) {
        unsigned int nMatches = player.dnaMatch(c, favors[c]);
        ratio += static_cast<double>(nMatches) / maxMatches[c];
    }
    ratio /= Life::NChromo;
    return ratio;
}

}
