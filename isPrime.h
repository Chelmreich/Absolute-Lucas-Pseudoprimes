#ifndef isPrime
#define isPrime

#include <inttypes.h>

typedef uint64_t uint64;

class Prime
{
        public:
                bool isPrime(uint64_t n);
                uint64_t mulMod(uint64_t a, uint64_t b, uint64_t m);
                uint64_t powMod(uint64_t a, uint64_t b, uint64_t m);
        private:
                bool sprp(uint64_t a, uint64_t b);
};
extern Prime isPrime;

#endif