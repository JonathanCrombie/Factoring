/* Public Domain.  See the LICENSE file.                            */

/* Prints a range of prime numbers                                  */

/* Super simple algorithm -- just runs the Sieve of Eratosthenes    */
/* for the full range from 0 to the last number of the range.       */

/* No dependencies.                                                 */
/* On linux, try:  cc prime_range.c -o prime_range                  */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

unsigned int  mask[32] = {0x00000001,0x00000002,0x00000004,0x00000008,
                          0x00000010,0x00000020,0x00000040,0x00000080,
                          0x00000100,0x00000200,0x00000400,0x00000800,
                          0x00001000,0x00002000,0x00004000,0x00008000,
                          0x00010000,0x00020000,0x00040000,0x00080000,
                          0x00100000,0x00200000,0x00400000,0x00800000,
                          0x01000000,0x02000000,0x04000000,0x08000000,
                          0x10000000,0x20000000,0x40000000,0x80000000 };


int64_t isqrt( int64_t number );

int main( int argc, char * argv[] ) {

  if ( argc != 3 ) {
    fprintf( stderr, "Usage: prime_range start end\n");
    return 1;
  }

  int64_t begin = atol( argv[1] );
  int64_t max_limit = 1000000000000000000L;
  if ( begin < 0 || begin > max_limit ) {
    fprintf( stderr, "Error: begin range must >= 0 and <= %ld. Aborting.\n\n", max_limit );
    return 1;
  }

  int64_t end = atol( argv[2] );
  if ( end < 0 || end > max_limit ) {
    fprintf( stderr, "Error: end range must >= 0 and <= %ld. Aborting.\n\n", max_limit );
    return 1;
  }

  if ( begin > end ) {
    fprintf( stderr, "Error: Begin range must be less than end range. Aborting.\n\n" );
    return 1;
  }

  // each bit in "array" corresponds to one whole number from 0 to limit
  // eg. the number  0 --> bit  0 --> array[0] & mask[0]
  //     the number 31 --> bit 31 --> array[0] & mask[31]
  //     the number 32 --> bit 32 --> array[1] & mask[0]
  //
  // after all primes are computed, a bit set to 1 will represent a non-prime
  // and a bit set to 0 will represent a prime

  uint32_t* array = (uint32_t *) calloc( end / 32 + 1, sizeof(uint32_t) );
  if ( array == NULL ) {
    fprintf( stderr, "Error: Failed to allocate memory. Aborting.\n\n" );
    return 1;
  }

  array[0] |= mask[0]; // corresponds to the number 0 which is not prime
  array[0] |= mask[1]; // corresponds to the number 1 which is not prime

  long sqrroot_of_upper_limit = isqrt(end) + 1;
  int64_t quot = 0;
  int64_t rem = 0;
  int64_t i = 0;
  int64_t j = 0;
  for (i = 2; i <= sqrroot_of_upper_limit; i++) {
    quot = i >> 5;
    rem = i & 0x0000001F;
    if (!(mask[rem] & array[quot])) {
      for (j = i+i; j <= end; j += i) {
        quot = j >> 5;
        rem = j & 0x0000001F;
        array[quot] |= mask[rem];
      }
    }
  }

  // Print out the primes
  for (i = begin; i <= end; i++) {
    quot = i >> 5;
    rem = i & 0x0000001F;
    if (!(mask[rem] & array[quot]))
      printf( "%ld\n", i );
  }

  if ( array != NULL ) {
    free( array );
    array = NULL;
  }

  return 0;
}

// Simple integer square root algorithm from Google AI search
int64_t isqrt( int64_t number ) {
  int64_t a = number;
  int64_t b = (number + 1) / 2; // Initial guess

  while (a > b) {
    a = b;
    b = (b + number / b) / 2;
  }

  // Ensure the result is the floor of the square root
  if (a * a > number)
    a--;

  return a;
}

