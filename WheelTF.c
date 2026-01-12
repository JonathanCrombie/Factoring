/* Public Domain.  See the LICENSE file. */


/* A simple program to illustrate trial factoring a number using the wheel   */
/* factorization method.                                                     */

/* To compile, the GMP library needs to be already installed.                */
/* See https://gmplib.org                                                    */
/* On linux, try:  cc WheelTF.c -lgmp -o WheelTF                             */


/* The maximum integer we attempt to trial factor is hard-coded to approx.   */
/* 4 x 10^9.  Note, this program is still very inefficient even for finding  */
/* factors of this size.  Not much attempt has been made at optimization.    */
/* For more "serious" factoring programs, check out for example msieve,      */
/* yafu, and gmp-ecm.                                                        */

/* Some test numbers:                                                         */
/* WheelTF 1022117                   --> 1009.1013                            */
/* WheelTF 100160063                 --> 10007.10009                          */
/* WheelTF 10002200057               --> 100003.100019                        */
/* WheelTF 1000036000099             --> 1000003.1000033                      */
/* WheelTF 100000980001501           --> 10000019.10000079                    */
/* WheelTF 10000004400000259         --> 100000007.100000037              ~1s */
/* WheelTF 1000000016000000063       --> 1000000007.1000000009            ~5s */



#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>

struct factor_info {
mpz_t           the_factor;
long            occurrences;
char            factor_status;
};

struct factor_infos {
long                 count;
struct factor_info*  the_factors;
};

char quickprimecheck( mpz_t );
void TFDivideOut( mpz_t, long, char*, mpz_t, struct factor_infos* );
void WheelTF( mpz_t, struct factor_infos* );
long ComputeOccurrences( mpz_t, mpz_t );
void Init_Factor_Infos( struct factor_infos* );
void AddFactorInfo( struct factor_infos*, mpz_t, long, char );
void Cleanup_Factor_Infos( struct factor_infos* );

// Explanation of where the numbers come from.

// eg. looking at wheel5, wheel5 means wheel of primes 2 ... 5 == 2,3,5
// wheel circumference = 2.3.5 = 30
// wheel spikes removed if not coprime with circumference.
// (ie. sieve of eratosthenes for 2 ... 5 up to the number 30)
// leaves spikes at 1,7,11,13,17,19,23,29
// differences are 6,4,2,4,2,4,6,2
// rotated differences to match starting at next prime == 7,
// yields 4,2,4,2,4,6,2,6

// wheel3 and wheel5 are just here for illustration purposes
const uint8_t wheel3[2] = { 2, 4 };
const uint8_t wheel5[8] = { 4, 2, 4, 2, 4, 6, 2, 6 };

// wheel7 is what we actually use
const uint8_t wheel7[48] = { 2, 4, 2, 4, 6, 2, 6, 4, 2, 4,
                             6, 6, 2, 6, 4, 2, 6, 4, 6, 8,
                             4, 2, 4, 2, 4, 8, 6, 4, 6, 2,
                             4, 6, 2, 6, 6, 4, 2, 4, 6, 2,
                             6, 4, 2, 4, 2, 10, 2, 10 };

int main( int argc, char * argv[] ) {

  if ( argc != 2 ) {
    printf( "\nUsage: WheelTF n\n\n" );
    return 1;
  }

  mpz_t n;
  mpz_init_set_str( n,  argv[1], 10 );

  if ( mpz_cmp_ui( n, 2 ) < 0 ) {
    printf("\nThe number must be >= 2.  Aborting.\n\n");
    mpz_clear(n);
    return 1;
  }

  struct factor_infos Factor_Infos;
  Init_Factor_Infos( &Factor_Infos );

  WheelTF( n, &Factor_Infos );

  printf("\n");
  long i = 0;
  for ( ; i < Factor_Infos.count; i++ ) {
    if ( i > 0 )
      printf( ".");

    gmp_printf( "%Zd", Factor_Infos.the_factors[i].the_factor );

    if ( Factor_Infos.the_factors[i].factor_status == 'C' )
      printf("C");

    if ( Factor_Infos.the_factors[i].occurrences > 1 )
      printf( "^%ld", Factor_Infos.the_factors[i].occurrences );
  }
  printf("\n");

  Cleanup_Factor_Infos( &Factor_Infos );
  mpz_clear( n );

  return 0;
}

// Returns number type.  'C' --> Composite, 'P' --> Prime (very probably), 'N' --> Neither (ie. the number 1)
char quickprimecheck( mpz_t the_number ) {

  char retval = 'C';
  if ( mpz_cmp_ui( the_number, 1 ) == 0 ) {
    retval = 'N';
    return retval;
  }

  int factor_type = mpz_probab_prime_p( the_number, 30 );
  if ( factor_type != 0 )
    retval = 'P';

  return retval;
}


// divide out denominator from running_N
void TFDivideOut( mpz_t running_N, long ldenominator, char* running_N_status, mpz_t square_root, struct factor_infos* Factor_Infos ) {

  mpz_t denominator;
  mpz_init_set_ui( denominator, ldenominator );

  long occurrences = ComputeOccurrences( running_N, denominator );

  long i = 0;
  for ( i = occurrences; i > 0; i-- )
    mpz_divexact_ui( running_N, running_N, ldenominator );

  *running_N_status = quickprimecheck( running_N );
  if ( *running_N_status == 'N' )
    mpz_set_ui( square_root, 1 );
  else
    mpz_sqrt( square_root, running_N );

  AddFactorInfo( Factor_Infos, denominator, occurrences, 'P' );

  mpz_clear( denominator );
}

// Compute the prime factorization of a general number
void WheelTF( mpz_t the_number, struct factor_infos* Factor_Infos ) {
  if ( Factor_Infos == NULL )
    return;

  // Wheel Factorization. see eg. http://programmingpraxis.com/2009/05/08/wheel-factorization/
  // Use a 2,3,5,7 wheel (or 'wheel7' for short)

  mpz_t running_N;
  mpz_init_set( running_N, the_number );
  char running_N_status = quickprimecheck( running_N );

  mpz_t square_root;
  mpz_init( square_root );
  mpz_sqrt( square_root, running_N );

  // This wheel only starts to be applied at number 11, so 2, 3, 5, and 7
  // will have to be manually checked first.
  if ( mpz_divisible_ui_p( running_N, 2 ) != 0 )
    TFDivideOut( running_N, 2, &running_N_status, square_root, Factor_Infos );

  if ( mpz_divisible_ui_p( running_N, 3 ) != 0 )
    TFDivideOut( running_N, 3, &running_N_status, square_root, Factor_Infos );

  if ( mpz_divisible_ui_p( running_N, 5 ) != 0 )
    TFDivideOut( running_N, 5, &running_N_status, square_root, Factor_Infos );

  if ( mpz_divisible_ui_p( running_N, 7 ) != 0 )
    TFDivideOut( running_N, 7, &running_N_status, square_root, Factor_Infos );

  // for speed reasons, we will use an unsigned long as the trial factor and as the tf upper limit.
  unsigned long tf_upperlimit =  4000000000ul;

  if ( mpz_cmp_ui( square_root, tf_upperlimit ) < 0 )
    tf_upperlimit = mpz_get_ui( square_root );

  uint8_t i = 0;
  unsigned long tf = 11;

  for ( i = 0; tf <= tf_upperlimit; tf += wheel7[i], i = ( i == 47 ? 0 : i + 1 ) ) {
      if ( mpz_divisible_ui_p( running_N, tf ) != 0 ) {
        TFDivideOut( running_N, tf, &running_N_status, square_root, Factor_Infos );
        if ( running_N_status != 'C' )
          break;

        if ( mpz_cmp_ui( square_root, tf_upperlimit ) < 0 )
          tf_upperlimit = mpz_get_ui( square_root );
      }
  }

  if ( mpz_cmp_ui( running_N, 1 ) != 0 )
    AddFactorInfo( Factor_Infos, running_N, 1, running_N_status );

  mpz_clear( square_root );
  mpz_clear( running_N );
}

// Compute the number of times denominator divides evenly into numerator
long ComputeOccurrences( mpz_t numerator, mpz_t denominator ) {

  // not handled values
  if ( mpz_cmp_ui( denominator, 1 ) <= 0 )
    return -1;

  long occurrences = 0;

  mpz_t quotient;
  mpz_init_set( quotient, numerator );

  while ( mpz_divisible_p( quotient, denominator ) ) {
    mpz_divexact( quotient, quotient, denominator );
    occurrences++;
  }

  mpz_clear( quotient );

  return occurrences;
}

// Initialize factor_infos
void Init_Factor_Infos( struct factor_infos* Factor_Infos ) {
  if ( Factor_Infos == NULL )
    return;
  Factor_Infos->count = 0;
  Factor_Infos->the_factors = NULL;
}

// Add a factor to factor_infos structure
void AddFactorInfo( struct factor_infos* Factor_Infos, mpz_t the_factor, long occurrences, char factor_status ) {

  long index = 0;

  // allocate memory
  if ( Factor_Infos->count == 0 ) {
    Factor_Infos->count = 1;
    Factor_Infos->the_factors = (struct factor_info*) calloc( 1, sizeof(struct factor_info) );
  }
  else {
    Factor_Infos->count++;
    Factor_Infos->the_factors = (struct factor_info*) realloc( Factor_Infos->the_factors, sizeof(struct factor_info) * Factor_Infos->count );
    index = Factor_Infos->count - 1;
    memset( &Factor_Infos->the_factors[index], 0, sizeof(struct factor_info) );
  }

  mpz_init( Factor_Infos->the_factors[index].the_factor );
  mpz_set( Factor_Infos->the_factors[index].the_factor, the_factor );
  Factor_Infos->the_factors[index].occurrences    = occurrences;
  Factor_Infos->the_factors[index].factor_status  = factor_status;
}

// Free the memory allocated
void Cleanup_Factor_Infos( struct factor_infos* Factor_Infos ) {
  if ( Factor_Infos == NULL )
    return;

  long i;
  for ( i = 0; i < Factor_Infos->count; i++ )
    mpz_clear( Factor_Infos->the_factors[i].the_factor );

  if ( Factor_Infos->the_factors != NULL ) {
    free( Factor_Infos->the_factors );
    Factor_Infos->the_factors = NULL;
  }
  Factor_Infos->count = 0;
}



