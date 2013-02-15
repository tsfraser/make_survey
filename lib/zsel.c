#pragma once
#ifndef ZSEL_INCLUDED
#define ZSEL_INCLUDED

#include <math.h>

#include <check_alloc.c>
#include <simple_reader.c>
#include <spline_gsl.c>

typedef struct {
    SPLINE *spl;
    double zmin;
    double zmax;
} ZSEL;

ZSEL *
zsel_init( void )
{
    ZSEL *zs;
    zs = ( ZSEL * ) check_alloc( 1, sizeof( ZSEL ) );

    zs->spl = spline_init( SPLINE_TYPE_AKIMA ); /* allow discontinuities, with potential deviations */
    zs->zmin = 1e100;
    zs->zmax = -1e100;

    return zs;
}

ZSEL *
zsel_kill( ZSEL * zs )
{

    spline_kill( zs->spl );
    CHECK_FREE( zs );
    return zs;
}

ZSEL *
zsel_read_file( char const *const filename )
{
    ZSEL *zs;
    double z, frac;
    simple_reader *sr;

    zs = zsel_init(  );

    sr = sr_init( filename );
    while( sr_readline( sr ) ) {
        char *line;
        int check;

        line = sr_line( sr );

        if( sr_line_isempty( sr ) )
            continue;
        if( line[0] == '#' )
            continue;

        check = sscanf( line, "%lf %lf", &z, &frac );
        if( check != 2 ) {
            fprintf( stderr,
                     "ZSEL Error: could not read two columns on line %d in file: %s\n",
                     sr_linenum( sr ), sr_filename( sr ) );
            exit( EXIT_FAILURE );
        }

        if( frac > 1.0 || frac < 0.0 ) {
            fprintf( stderr,
                     "ZSEL Error: Input requires 0.0 <= f <= 1.0 (not f = %g) on line %d in file: %s\n",
                     frac, sr_linenum( sr ), sr_filename( sr ) );
            exit( EXIT_FAILURE );
        }

        /* find zmin,zmax */
        zs->zmin = fmin( zs->zmin, z );
        zs->zmax = fmax( zs->zmax, z );

        spline_data_add( zs->spl, z, frac );
    }
    sr_kill( sr );

    return zs;
}

static inline double
zsel_eval( ZSEL const *const zs, const double z )
{
    /* gracefully handle boundaries, then continue */
    if( z < zs->zmin )
        return 0.0;
    if( z > zs->zmax )
        return 0.0;

    return spline_eval( zs->spl, z );
}

#endif
