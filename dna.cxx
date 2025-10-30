/*
 * Sin-Yaw Wang <swang24@scu.edu>
 */
#include <bitset>
#include "rng.h"
#include "dna.h"

namespace csen79 {

// Student to implement all of these member functions
DNA::DNA(unsigned int v) {
    codes = static_cast<Gene>(v) & static_cast<Gene>(mask);
}

// Read bit i (0 = least significant). Throws if i is out of range.
bool DNA::getCode(const int i) const {
    if (i < 0 || i >= nCode) throw std::out_of_range("DNA::getCode index out of range");
    return (static_cast<unsigned int>(codes) >> i) & 0x1u;
}

// Write bit i to v (true->1, false->0). Throws if i is out of range.
void DNA::setCode(const int i, bool v) {
    if (i < 0 || i >= nCode) throw std::out_of_range("DNA::setCode index out of range");
    const unsigned int bit = 1u << i;
    unsigned int tmp = static_cast<unsigned int>(codes);
    if (v) tmp |= bit;
    else   tmp &= ~bit;
    codes = static_cast<Gene>(tmp & mask);
}

// Constant size (number of genetic codes).
size_t DNA::size() const { return static_cast<size_t>(nCode); }

// Number of 1-bits currently set.
size_t DNA::count() const {
#if defined(__cpp_lib_bitops) || (defined(__GNUC__) && !defined(__clang__))
    // Use builtin popcount when available
    return static_cast<size_t>(__builtin_popcount(static_cast<unsigned int>(codes)));
#else
    // Portable approach via std::bitset
    return std::bitset<nCode>(static_cast<unsigned long long>(codes)).count();
#endif
}

// How many favored bits match between two DNAs.
// We define "match" as positions where both have 1 (intersection), so
// matchDNA = popcount(this->codes & other.codes).
unsigned int DNA::matchDNA(const DNA& other) const {
#if defined(__cpp_lib_bitops) || (defined(__GNUC__) && !defined(__clang__))
    return static_cast<unsigned int>(__builtin_popcount(
        static_cast<unsigned int>(codes & other.codes)
    ));
#else
    return static_cast<unsigned int>(
        std::bitset<nCode>(static_cast<unsigned long long>(codes & other.codes)).count()
    );
#endif
}

}   // namespace csen79
