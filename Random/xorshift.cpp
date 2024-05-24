/*
* Updates the RNG state in cooperation with in - warp neighbors .
* Uses a block of shared memory of size
* ( WARPSIZE + WORDSHIFT + 1) * NWARPS + WORDSHIFT + 1.
* Parameters :
* state : RNG state
* tid: thread index in block
* stateblock : shared memory block for states
* Returns :
 * updated state
 */
 __device__ state_t rng update( state_t state , int tid ,
 v o l a t i l e state_t * stateblock )
 {

 /* Indices . */
 int wid = tid / WARPSIZE ; // Warp index in block
 int lid = tid % WARPSIZE ; // Thread index in warp
 int woff = wid * ( WARPSIZE + WORDSHIFT + 1) + WORDSHIFT + 1;
 // warp offset
 /* Shifted indices . */
 int lp = lid + WORDSHIFT ; // Left word shift
 int lm = lid - WORDSHIFT ; // Right word shift

 /* << A. */
 stateblock [ woff + lid ] = state ; // Share states
 state ^= stateblock [ woff + lp ] << RAND_A ; // Left part
 state ^= stateblock [ woff + lp + 1] >> WORD - RAND_A ; // Right part

 /* >> B. */
 stateblock [ woff + lid ] = state ; // Share states
 state ^= stateblock [ woff + lm - 1] << WORD - RAND_B ; // Left part
 state ^= stateblock [ woff + lm ] >> RAND_B ; // Right part

 /* << C. */
 stateblock [ woff + lid ] = state ; // Share states
 state ^= stateblock [ woff + lp ] << RAND_C ; // Left part
 state ^= stateblock [ woff + lp + 1] >> WORD - RAND_C ; // Right part

 return state ;
 }

