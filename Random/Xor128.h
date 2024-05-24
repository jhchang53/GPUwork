#ifndef _inc_Xor128
#define _inc_Xor128

/*  this file is needed to call func. of Xor128.cpp  */
#include <stdint.h>
class Xor128 {
public:
  Xor128();
  ~Xor128();
  double ranf();
  uint64_t next_int();
  void jump(void);
  uint64_t next(void);
  /*  isotropic angle  */
  void direcS2(double w[]);
private:
  
  uint64_t seed_state;
  uint64_t rand_s[2];
};
#endif
