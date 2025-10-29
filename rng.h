/*
 * Sin-Yaw Wang <swang24@scu.edu>
 *
 * Simple random number generation helpers.
 * Implementation in main.cxx
 */
#ifndef RNG_H
#define RNG_H

// Returns random integer in [lo, hi]
int rand_int(int lo, int hi);

// Returns true/false with 50% probability
bool rand_bool();

#endif // RNG_H
