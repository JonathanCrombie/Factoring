/* Public Domain.  See the LICENSE file.                                     */

/* The simplest possible implementation of Pollard's Rho Factorization       */
/* algorithm.                                                                */
/* https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm                   */

/* To compile, the GMP library needs to be already installed.                */
/* See https://gmplib.org                                                    */
/* On linux, try:  g++ rho.cpp -lgmp -o rho                                  */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void g( mpz_t, mpz_t );

int main( int argc , char * argv[] ) {

  if ( argc != 2 ) {
    printf("\nUsage: rho N\n");
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

  if ( mpz_probab_prime_p( n, 30 ) ) {
    printf( "N is a probable prime. Aborting.\n\n" );
    mpz_clear( n );
    return 1;
  }

  mpz_t x;
  mpz_init( x );
  mpz_t y;
  mpz_init( y );
  mpz_t d;
  mpz_init( d );

  mpz_set_ui( d, 1 );
  mpz_set_ui( x, 2 );
  mpz_set( y, x );


  mpz_t tempZ1;
  mpz_init( tempZ1 );

  while ( !mpz_cmp_ui( d, 1 ) ) {
    g( x, n );
    g( y, n );
    g( y, n );

    mpz_sub( tempZ1, x, y );
    mpz_abs( tempZ1, tempZ1 );
    mpz_gcd( d, tempZ1, n );
  }

  if ( !mpz_cmp( d, n ) )
    printf( "Failure\n" );
  else
    gmp_printf( "Found a non-trivial factor: %Zd\n", d );

  mpz_clear( tempZ1 );
  mpz_clear( d );
  mpz_clear( y );
  mpz_clear( x );

  return 0;
 }

// computes the polynomial in x which will be hard coded to "x^2 + 1" mod n
void g( mpz_t x, mpz_t n ) {
  mpz_mul( x, x, x );
  mpz_add_ui( x, x, 1 );
  mpz_mod( x, x, n );
}

