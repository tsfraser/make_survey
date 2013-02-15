/* Simplified wrapper around GSL functionality */

#pragma once
#ifndef SPLINE_INCLUDED
#define SPLINE_INCLUDED

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#include <check_alloc.c>
#include <simple_array.c>
#include <assert.h>

#define SPLINE_TYPE_LINEAR gsl_interp_linear
#define SPLINE_TYPE_CUBIC gsl_interp_cspline
#define SPLINE_TYPE_AKIMA gsl_interp_akima

typedef struct {
    gsl_interp *interp;
    gsl_interp_type const *type;
    gsl_interp_accel *accel;
    size_t len;
    simple_array *xsa;
    simple_array *ysa;
} SPLINE;

static inline int
spline_check_not_null( void const *const v )
{
    if( NULL == v ) {
        fprintf( stderr, "SPLINE Error: expected NON NULL value!\n" );
        assert( NULL == v );
    }
    return TRUE;
}

void
spline_alloc( SPLINE * s, size_t data_len_guess )
{
    spline_check_not_null( s );
    s->len = 0;
    s->interp = NULL;
    s->accel = NULL;
    s->xsa = sa_init( data_len_guess, sizeof( double ) );
    s->ysa = sa_init( data_len_guess, sizeof( double ) );
}

SPLINE *
spline_init( void const *type )
{
    SPLINE *s;

    s = ( SPLINE * ) check_alloc( 1, sizeof( SPLINE ) );
    s->type = ( gsl_interp_type const * ) type;
    spline_alloc( s, 1000 );

    return s;
}

void
spline_clean( SPLINE * s )
{
    spline_check_not_null( s );
    gsl_interp_free( s->interp );
    gsl_interp_accel_free( s->accel );
    s->interp = NULL;
    s->accel = NULL;

    s->xsa = sa_kill( s->xsa );
    s->ysa = sa_kill( s->ysa );
}

SPLINE *
spline_kill( SPLINE * s )
{
    spline_clean( s );
    CHECK_FREE( s );
    return NULL;
}

static inline void
spline_data_add( SPLINE * s, const double x, const double y )
{
    size_t len;
    double *xa, *ya;

    spline_check_not_null( s );
    len = s->len;

    if( s->interp != NULL ) {
        fprintf( stderr, "SPLINE Error: cannot add data after finalizing!\n" );
        exit( EXIT_FAILURE );
    }

    xa = ( double * ) sa_ensure_length( s->xsa, len + 1 );
    ya = ( double * ) sa_ensure_length( s->ysa, len + 1 );

    /* must be strictly ordered in increasing x */
    if( len > 0 )
        if( xa[len - 1] > x ) {
            fprintf( stderr, "SPLINE Error: input data must be in increasing value of x!\n" );
            exit( EXIT_FAILURE );
        }

    xa[len] = x;
    ya[len] = y;
    s->len += 1;
}

static inline void
spline_data_finalize( SPLINE * const s )
{
    double *xa, *ya;

    spline_check_not_null( s );
    /* first freeze X/Y arrays */
    xa = ( double * ) sa_set_length( s->xsa, s->len );
    ya = ( double * ) sa_set_length( s->ysa, s->len );

    /* actually initialize spline now */
    s->interp = gsl_interp_alloc( s->type, s->len );
    gsl_interp_init( s->interp, xa, ya, s->len );
    s->accel = gsl_interp_accel_alloc(  );
}

static inline double
spline_eval( SPLINE * s, const double xv )
{
    double *xa, *ya;

    spline_check_not_null( s );

    if( NULL == s->interp )
        spline_data_finalize( s );

    xa = ( double * ) sa_data( s->xsa );
    ya = ( double * ) sa_data( s->ysa );

    return gsl_interp_eval( s->interp, xa, ya, xv, s->accel );
}

static inline double
spline_deriv( SPLINE * const s, const double xv )
{
    double *xa, *ya;

    spline_check_not_null( s );

    if( NULL == s->interp )
        spline_data_finalize( s );

    xa = ( double * ) sa_data( s->xsa );
    ya = ( double * ) sa_data( s->ysa );

    return gsl_interp_eval_deriv( s->interp, xa, ya, xv, s->accel );
}

static inline double
spline_deriv2( SPLINE * const s, const double xv )
{
    double *xa, *ya;

    spline_check_not_null( s );

    if( NULL == s->interp )
        spline_data_finalize( s );

    xa = ( double * ) sa_data( s->xsa );
    ya = ( double * ) sa_data( s->ysa );

    return gsl_interp_eval_deriv2( s->interp, xa, ya, xv, s->accel );
}

static inline double
spline_integ( SPLINE * const s, const double a, const double b )
{
    double *xa, *ya;

    spline_check_not_null( s );

    if( NULL == s->interp )
        spline_data_finalize( s );

    xa = ( double * ) sa_data( s->xsa );
    ya = ( double * ) sa_data( s->ysa );

    return gsl_interp_eval_integ( s->interp, xa, ya, a, b, s->accel );
}

#endif
