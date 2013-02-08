/* Simplified wrapper that separates out using GSL RNG */

#pragma once
#ifndef RNG_INCLUDED
#define RNG_INCLUDED

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* initialize GSL_RNG and return a variable which holds state */

static inline void *
rng_init( const unsigned long int seed )
{
    gsl_rng *rng;

    /* this sets the type, MT19937 is a good RNG */
    rng = gsl_rng_alloc( gsl_rng_mt19937 );
    gsl_rng_set( rng, seed );

    return ( void * ) rng;
}

/* clean up to deallocate RNG structure */
static inline void *
rng_kill( void const *state )
{
    gsl_rng_free( ( gsl_rng * ) state );
    return NULL;
}

/* a wrapper to yield a uniform deviate */
static inline double
rng_uniform( void const *state )
{
    return gsl_rng_uniform( ( gsl_rng * ) state );
}

/* a wrapper to yield a gaussian deviate */
static inline double
rng_gaussian( void const *state, const double sigma )
{
    return gsl_ran_gaussian( ( gsl_rng * ) state, sigma );
}

#endif
