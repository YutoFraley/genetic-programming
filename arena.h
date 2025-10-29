/*
 * Sin-Yaw Wang <swang24@scu.edu>
 */
#ifndef ARENA_H
#define ARENA_H

namespace csen79 {

class Arena {
    static constexpr int NPop{1000};    // Arbitrary
    static constexpr int EVENTPERROUND{10};// Arbitrary
    std::vector<Life> population;
    Select competition;
public:
    static constexpr double MATURITY{0.9};
    Arena() {
        population.reserve(NPop);
        for (int i = 0; i < NPop; i++)
            population.push_back(Life());
    }
    void compete();
    void evolve();
    double matured();
    int highGeneration() {
        int highGen{0};
        for (auto& individual : population)
            if (individual.getGeneration() > highGen)
                highGen = individual.getGeneration();
        return highGen;
    }
    double averageScore() {
        long total{0};
        for (auto& individual : population)
            total += individual.getScore();
        return static_cast<double>(total) / NPop;
    }
};

}   // namespace csen79


#endif // ARENA_H