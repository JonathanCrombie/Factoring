/* Public Domain.  See the LICENSE file.                            */

/* Sieve of Eratosthenes                                            */

/* https://wikipedia.org/wiki/Sieve_of_Eratosthenes                 */
/* https://t5k.org/howmany.html#table (For prime counts)            */

/* No dependencies.                                                 */
/* On linux, try:  cc eratosthenes.c -o eratosthenes                */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

  printf( "\n" );
  printf( "         Sieve Of Eratosthenes\n" );
  printf( "\n" );
  printf( "\n" );
  printf( "\n" );
  printf( "Usage: eratosthenes limit\n" );
  printf( "\n" );
  printf( "\n" );
  printf( "NOTE: Memory usage in bytes will be limit / 8.\n" );
  printf( "\n" );
  printf( "eg. \"./eratosthenes 1000000000\" will use 125000000 bytes or apprx 120 MiB\n" );
  printf( "\n" );
  printf( "\n" );

  if ( argc != 2 ) {
    fprintf( stderr, "Usage: eratosthenes limit\n");
    return 1;
  }

  int64_t limit = atol( argv[1] );
  int64_t max_limit = 1000000000000000000L;
  if ( limit < 2 || limit > max_limit  ) {
    fprintf( stderr, "Error: limit must >= 2 and <= %ld. Aborting.\n\n", max_limit );
    return 1;
  }

  struct timespec  time_t0;
  clock_gettime(CLOCK_REALTIME, &time_t0);


  // each bit in "array" corresponds to one whole number from 0 to limit
  // eg. the number  0 --> bit  0 --> array[0] & mask[0]
  //     the number 31 --> bit 31 --> array[0] & mask[31]
  //     the number 32 --> bit 32 --> array[1] & mask[0]
  //
  // after all primes are computed, a bit set to 1 will represent a non-prime
  // and a bit set to 0 will represent a prime

  uint32_t* array = (uint32_t *) calloc( limit / 32 + 1, sizeof(uint32_t) );
  if ( array == NULL ) {
    fprintf( stderr, "Error: Failed to allocate memory. Aborting.\n\n" );
    return 1;
  }

  struct timespec  time_t1;
  clock_gettime(CLOCK_REALTIME, &time_t1);

  array[0] |= mask[0]; // corresponds to the number 0 which is not prime
  array[0] |= mask[1]; // corresponds to the number 1 which is not prime

  int64_t j;
  long sqrroot_of_limit = isqrt(limit) + 1;
  int64_t i = 2;
  int64_t quot = 0;
  int64_t rem = 0;
  for (; i <= sqrroot_of_limit; i++) {
    quot = i >> 5;  // dividing by 2^5 which is 32
    rem = i & 0x0000001F; // last 5 bits are the remainder after dividing by 32
    if (!(mask[rem] & array[quot])) {
      for (j = i+i; j <= limit; j += i) {
        quot = j >> 5; // dividing by 2^5 which is 32
        rem = j & 0x0000001F; // last 5 bits are the remainder after dividing by 32
        array[quot] |= mask[rem];
      }
    }
  }

  struct timespec  time_t2;
  clock_gettime(CLOCK_REALTIME, &time_t2);

  intmax_t elapsed_secs = time_t1.tv_sec - time_t0.tv_sec;
  intmax_t elapsed_nsecs = time_t1.tv_nsec - time_t0.tv_nsec;

  if ( elapsed_nsecs < 0 ) {
    elapsed_nsecs += 1000000000;
    --elapsed_secs;
  }
  printf( "\n" );
  printf( "Time To allocate memory (secs):   %jd.%09jd\n", elapsed_secs, elapsed_nsecs );


  elapsed_secs = time_t2.tv_sec - time_t1.tv_sec;
  elapsed_nsecs = time_t2.tv_nsec - time_t1.tv_nsec;

  if ( elapsed_nsecs < 0 ) {
    elapsed_nsecs += 1000000000;
    --elapsed_secs;
  }
  printf( "\n" );
  printf( "Time To compute primes  (secs):   %jd.%09jd\n", elapsed_secs, elapsed_nsecs );

  int64_t count = 0;
  for (i = 0; i <= limit; i++) {
    quot = i >> 5; // dividing by 2^5 which is 32
    rem = i & 0x0000001F; // last 5 bits are the remainder after dividing by 32
    if (!(mask[rem] & array[quot]))
      count++;
  }

  struct timespec  time_t3;
  clock_gettime(CLOCK_REALTIME, &time_t3);

  printf( "\n" );
  printf( "Total number of primes generated: %ld\n", count );


  elapsed_secs = time_t3.tv_sec - time_t2.tv_sec;
  elapsed_nsecs = time_t3.tv_nsec - time_t2.tv_nsec;

  if ( elapsed_nsecs < 0 ) {
    elapsed_nsecs += 1000000000;
    --elapsed_secs;
  }

  printf( "\n" );
  printf( "Time To count primes    (secs):   %jd.%09jd\n", elapsed_secs, elapsed_nsecs );

  printf( "\n" );

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

