/* Public Domain.  See the LICENSE file.                            */

/* Prints a range of prime numbers                                  */

/* Just prime_range.c, but won't compute primes if not necessary    */

/* No dependencies.                                                 */
/* On linux, try:  cc prime_range2.c -o prime_range2                */


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
    fprintf( stderr, "Usage: prime_range2 start end\n");
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

  // 2 chunks of allocated memory. 
  // First chunk holds the always required calculated primes from 2 to square root of the biggest number.
  // The second chunk will skip primes that are greater than where the first chunk ended and less than
  // the "begin" value. As usual, it will end the "end" value.
  int64_t sqrroot_of_upper_limit = isqrt(end) + 1;
  int64_t chunk1_size = sqrroot_of_upper_limit / 32 + 1;
  uint32_t* chunk1 = (uint32_t *) calloc( chunk1_size, sizeof(uint32_t) );
  if ( chunk1 == NULL ) {
    fprintf( stderr, "Error: Failed to allocate memory (chunk1). Aborting.\n\n" );
    return 1;
  }
  chunk1[0] |= mask[0]; // corresponds to the number 0 which is not prime
  chunk1[0] |= mask[1]; // corresponds to the number 1 which is not prime
  int64_t chunk1_end = chunk1_size * 32;

  int64_t chunk2_begin = 0;
  int64_t chunk2_end = 0;
  uint32_t* chunk2 = NULL;
  int64_t chunk2_size = 0;
  int64_t gap_size = 0;
  if ( end > chunk1_end ) { // need to allocate
    if ( begin - chunk1_end >= 32 ) // there is a gap
      gap_size = (begin - chunk1_end) / 32;
    chunk2_begin = chunk1_end + gap_size * 32;
    chunk2_end = (end / 32) * 32 + 32;
    chunk2_size = (chunk2_end - chunk2_begin) / 32;
    chunk2 = (uint32_t *) calloc( chunk2_size, sizeof(uint32_t) );
    if ( chunk2 == NULL ) {
      fprintf( stderr, "Error: Failed to allocate memory (chunk2). Aborting.\n\n" );
      return 1;
    }
  }

  int64_t quot = 0;
  int64_t rem = 0;
  int64_t i = 0;
  int64_t j = 0;
  int64_t tempz1 = 0;
  int64_t chunk2_first_i = 0;
  int64_t chunk2_max_offset = chunk2_size * 32 + (end - chunk2_end);
  for (i = 2; i <= sqrroot_of_upper_limit; i++) {
    quot = i >> 5;
    rem = i & 0x0000001F;
    if (!(mask[rem] & chunk1[quot])) {
      for (j = i+i; j <= chunk1_end; j += i) {
        quot = j >> 5;
        rem = j & 0x0000001F;
        chunk1[quot] |= mask[rem];
      }

      if ( chunk2 != NULL ) {
        tempz1 = ((chunk2_begin / i) * i) - chunk2_begin; // calculate offset into chunk2
        if ( tempz1 < 0 )
            tempz1 += i;

        for ( j = tempz1; j <= chunk2_max_offset; j += i ) {
          quot = j >> 5;
          rem = j & 0x0000001F;
          chunk2[quot] |= mask[rem];
        }
      }
    }
  }

  // print chunk1
  if ( begin <= chunk1_end ) {
    int64_t print_begin = begin;
    int64_t print_end =  end <= chunk1_end ? end : chunk1_end;
    for (i = print_begin; i <= print_end; i++) {
      quot = i >> 5;
      rem = i & 0x0000001F;
      if (!(mask[rem] & chunk1[quot]))
        printf( "%ld\n", i );
    }
  }

  // print chunk2
  if ( chunk2 != NULL ) {
    i =  begin > chunk2_begin ? begin - chunk2_begin : 0;
    j = chunk2_begin + i;
    for (; i <= chunk2_max_offset; i++,j++) {
      quot = i >> 5;
      rem = i & 0x0000001F;
      if (!(mask[rem] & chunk2[quot])) {
        printf( "%ld\n", j );
      }
    }
  }

  if ( chunk2 != NULL ) {
    free( chunk2 );
    chunk2 = NULL;
  }

  if ( chunk1 != NULL ) {
    free( chunk1 );
    chunk1 = NULL;
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

