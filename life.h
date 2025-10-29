#ifndef LIFE_H
#define LIFE_H

namespace csen79 {
/*
 * This abstracts the lifeform that share the same DNA, or chromosome structure.
 * 
 * A derived class can change the number of chromosomes, but not the length of the DNA.
 * For that, it must derive a new DNA class first.
 */
class Life {
public:
    static constexpr int NChromo = 2;
    Life cross(const Life& other) const;    // beget an offspring
    void newLife() { score = 0; }
    void newScore(long s) { score += s; }
    long getScore() const { return score; }
    const int getGeneration() const { return generation; }
    const unsigned int dnaMatch(int c, const DNA& other) const;
private:
    DNA dna[NChromo];
    int generation{0};
    long score{0};
};

} // namespace csen79

#endif // LIFE_H
