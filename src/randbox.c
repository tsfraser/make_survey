#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <rng_gsl.c>

int
main( int argc, char *argv[] )
{
    /*---Mass-distribution-variables-------*/
    double lbox, density;
    double x[3];
    unsigned long int random_seed;
    unsigned int i, k, np;
    void *rng;

    if( argc < 4 ) {
        fprintf( stderr, "Usage: %s  LBOX  NUM_PART  RANDOM_SEED  >  OUT_XYZ\n", argv[0] );
        exit( EXIT_FAILURE );
    }

    sscanf( argv[1], "%lf", &lbox );
    sscanf( argv[2], "%u", &np );
    sscanf( argv[3], "%lu", &random_seed );

    density = ( double ) np / ( lbox * lbox * lbox );

    fprintf( stderr, "%s> creating random box: lbox = %g\n", argv[0], lbox );
    fprintf( stderr, "%s> %u particles ;  density = %e \n", argv[0], np, density );
    fprintf( stderr, "%s> initializing RNG using seed = %lu\n", argv[0], random_seed );

    rng = rng_init( random_seed );
    for( i = 0; i < np; i++ ) {
        for( k = 0; k < 3; k++ )
            x[k] = lbox * rng_uniform( rng );
        fprintf( stdout, "%12.6f %12.6f %12.6f\n", x[0], x[1], x[2] );
    }
    rng_kill( rng );

    return EXIT_SUCCESS;
}
