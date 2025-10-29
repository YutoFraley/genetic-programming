/*
 * Sin-Yaw Wang <swang24@scu.edu>
 */
#include <bitset>
#include <stdexcept>
#include <algorithm>
#include "rng.h"
#include "dna.h"
#include "life.h"

namespace csen79 {

// student to rewrite both
Life Life::cross(const Life& other) const {
    Life child;
    child.generation = std::max(this->generation, other.generation) + 1;
    
    // For each chromosome, randomly select genes from parents with possible mutation
    for (int c = 0; c < NChromo; c++) {
        // For each bit in the DNA
        for (size_t i = 0; i < child.dna[c].size(); i++) {
            // Small chance of mutation (1%)
            if (rand_int(1, 100) == 1) {
                child.dna[c].setCode(i, rand_bool());
            } else {
                // Randomly pick from one parent or the other
                if (rand_bool()) {
                    child.dna[c].setCode(i, this->dna[c].getCode(i));
                } else {
                    child.dna[c].setCode(i, other.dna[c].getCode(i));
                }
            }
        }
    }
    
    return child;
}

const unsigned int Life::dnaMatch(int c, const DNA& other) const {
    if (c < 0 || c >= NChromo)
        throw std::out_of_range("Life::dnaMatch: chromosome index out of range");
    return dna[c].matchDNA(other);
}

}   // namespace csen79
