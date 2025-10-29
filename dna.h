/*
 * Sin-Yaw Wang <swang24@scu.edu>
 */
#ifndef DNA_H
#define DNA_H

namespace csen79 {
class DNA {
    static constexpr int nCode = 16;
    static constexpr unsigned int mask = (1UL << nCode) - 1;
public:
    using Gene = unsigned short;    // are you in agreement with this choice?
    DNA(): DNA(rand_int(0, mask)) {}    // default is to have random DNA codes
    DNA(unsigned int v);
    bool getCode(const int i) const;
    void setCode(const int i, bool v);
    size_t size() const;    // length
    size_t count() const;   // count of 1s
    unsigned int matchDNA(const DNA& other) const;  // how many bits that match with the other
private:
    Gene codes;
};


} // namespace csen79

#endif // DNA_H
