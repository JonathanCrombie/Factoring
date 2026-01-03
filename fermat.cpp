/* Public Domain.  See the LICENSE file. */

/* A basic implementation of Fermat's Factoring Method straight from the     */
/* Wikipedia article.                                                        */
/* https://en.wikipedia.org/wiki/Fermat%27s_factorization_method             */

/* To compile, the GMP library needs to be already installed.                */
/* See https://gmplib.org                                                    */
/* On linux, try:  g++ fermat.cpp -lgmp -o fermat                            */

/* eg. try: ./fermat 5959                                                    */


#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>



int main(int argc , char * argv[]) {

  if (argc != 2) {
    printf("\nUsage: fermat N\n");
    return 1;
  }

  mpz_t n;
  mpz_init( n );

  mpz_set_str( n, argv[1], 10 );

  if ( mpz_cmp_ui( n, 100 ) < 0 ) {
    printf( "Lower bound on N is currently 100. Aborting.\n\n" );
    mpz_clear( n );
    return 1;
  }

  if ( mpz_divisible_2exp_p( n, 1 ) ) {
    printf( "N must be an odd number. Aborting.\n\n" );
    mpz_clear( n );
    return 1;
  }

  if ( mpz_probab_prime_p( n, 30 ) ) {
    printf( "N is a probable prime. Aborting.\n\n" );
    mpz_clear( n );
    return 1;
  }

  if ( mpz_perfect_square_p( n ) ) {
    printf( "N is a perfect square\n\n" );
    mpz_t sqrroot;
    mpz_init( sqrroot );
    mpz_sqrt( sqrroot, n );
    gmp_printf( "N = %Zd * %Zd\n\n", sqrroot, sqrroot );
    mpz_clear( sqrroot );
    mpz_clear( n );
    return 0;
  }

  mpz_t a;
  mpz_init( a );
  mpz_sqrt( a, n );
  mpz_add_ui( a, a, 1 );

  mpz_t b2; // b squared
  mpz_init( b2 );

  mpz_t tempZ1;
  mpz_init( tempZ1 );
  mpz_mul( tempZ1, a, a );
  mpz_sub( b2, tempZ1, n );

  while ( !mpz_perfect_square_p( b2 ) ) {
    mpz_mul_2exp( tempZ1, a, 1 );
    mpz_add( b2, b2, tempZ1 );
    mpz_add_ui( b2, b2, 1 );
    mpz_add_ui( a, a, 1 );
  }

  mpz_t b;
  mpz_init( b );
  mpz_sqrt( b, b2 );

  mpz_t a_minus_b;
  mpz_init( a_minus_b );
  mpz_t a_plus_b;
  mpz_init( a_plus_b );

  mpz_sub( a_minus_b, a, b );
  mpz_add( a_plus_b, a, b );

  gmp_printf( "%Zd %Zd\n", a_plus_b, a_minus_b );

  mpz_clear( a_plus_b );
  mpz_clear( a_minus_b );
  mpz_clear( b );
  mpz_clear( tempZ1 );
  mpz_clear( b2 );
  mpz_clear( a );
  mpz_clear( n );

  return 0;
 }
