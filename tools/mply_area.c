#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include <minimal_mangle.c>

int
main( int argc, char **argv )
{
    MANGLE_PLY *ply;
    double min_weight = 0.0;
    double area, warea;
    const double str2sdeg = ( 180.0 * 180.0 / M_PI / M_PI );

    if( argc < 2 ) {
        printf( "Usage: %s  POLYGON  [MIN_WEIGHT]\n", argv[0] );
        return EXIT_FAILURE;
    }

    fprintf( stderr, "Reading polygon file: %s\n", argv[1] );
    ply = mply_read_file( argv[1] );

    if( argc > 2 ) {
        min_weight = strtod( argv[2], NULL );
    }

    fprintf( stderr, "Filtering: weight >= %g\n", min_weight );

    area = mply_area_total( ply, min_weight );
    fprintf( stdout, "Total area: \n  str:   %g\n  deg^2: %g \n", area, area * str2sdeg );
    warea = mply_area_weighted_total( ply, min_weight );
    fprintf( stdout, "Weighted area: \n  str:   %g\n  deg^2: %g \n", warea, warea * str2sdeg );

    ply = mply_kill( ply );

    return EXIT_SUCCESS;
}
