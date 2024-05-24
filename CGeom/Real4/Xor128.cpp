/*  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */


/* This is xoroshiro128++ 1.0, one of our all-purpose, rock-solid,
   small-state generators. It is extremely (sub-ns) fast and it passes all
   tests we are aware of, but its state space is large enough only for
   mild parallelism.

   For generating just floating-point numbers, xoroshiro128+ is even
   faster (but it has a very mild bias, see notes in the comments).

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */
#include <stdint.h>
#include <math.h>
#include "Xor128.h"


static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

Xor128::Xor128()
{
  seed_state = 0;
  rand_s[0] = next_int();
  rand_s[1] = next_int();
};	// Xor128::Xor128

Xor128::~Xor128()
{
};

double Xor128::ranf()
{
   uint64_t ranu = next();
   return ranu*5.42101086242752217e-20;
};

void Xor128::direcS2(float w[])
{
  /*  uniform direction cosines of 3D sphere  */
  double v1,v2;
  double S = 99.9;
  while(S > 1.0) {
    v1 = 2*ranf()-1.0;
    v2 = 2*ranf()-1.0;
    S = v1*v1+v2*v2;
  }
  float rtS = sqrt(1.0-S);
  w[0] = 2*v1*rtS;
  w[1] = 2*v2*rtS;
  w[2] = 1.0-2*S;
};	// Xor128::direcS2


// uint64_t rand_s[2];

uint64_t Xor128::next(void) {
	const uint64_t s0 = rand_s[0];
	uint64_t s1 =rand_s[1];
	const uint64_t result = rotl(s0 + s1, 17) + s0;

	s1 ^= s0;
	rand_s[0] = rotl(s0, 49) ^ s1 ^ (s1 << 21); // a, b
	rand_s[1] = rotl(s1, 28); // c

	return result;
}


/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void Xor128::jump(void) {
	static const uint64_t JUMP[] = { 0x2bd7a6a6e99c2ddc, 0x0992ccaf6a6fca05 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= rand_s[0];
				s1 ^= rand_s[1];
			}
			next();
		}

	rand_s[0] = s0;
	rand_s[1] = s1;
}

/* Splitmix64 is the default pseudo-random number generator algorithm in Java and is included / available in many other languages */
 /* The state can be seeded with any (upto) 64 bit integer value. */
// uint64_t seed_state;

uint64_t  Xor128::next_int() {
  seed_state += 0x9e3779b97f4a7c15;  /* increment the state variable */
  uint64_t z = seed_state;		/* copy the state to a working variable */
  /* xor the variable with the variable right bit shifted 30 then multiply by a constant */
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  /* xor the variable with the variable right bit shifted 27 then multiply by a constant */
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);	/* return the variable xored with itself right bit shifted 31 */
};

