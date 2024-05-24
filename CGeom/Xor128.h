#ifndef _inc_Xor128
#define _inc_Xor128

#include <stdint.h>
extern uint64_t xor128_seed_state;
extern uint64_t xor128_s[2];

extern void init_xor128();
extern void jump_xor128(void);
extern double ranf();
#endif
