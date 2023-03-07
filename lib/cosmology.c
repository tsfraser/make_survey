#pragma once
#ifndef COSMOLOGY_INCLUDED
#define COSMOLOGY_INCLUDED

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_const_mksa.h>

#include <check_alloc.c>

#define COSMO_INT_SIZE 1000

/* we make assumptions on units below: Mpc and km/s */
#define SPEED_OF_LIGHT (GSL_CONST_MKSA_SPEED_OF_LIGHT / 1000.0)
#define HUBBLE_UNITS 100.0

typedef struct {
    double omega_m;
    double omega_l;
    double h;
    double w0;
    double wa;
    gsl_integration_workspace *work;
    gsl_function *df;
} COSMO;

static inline void
cosmo_alloc( COSMO * const c )
{
    /* initialize GSL integration routines */
    c->work = gsl_integration_workspace_alloc( COSMO_INT_SIZE );
    c->df = ( gsl_function * ) check_alloc( 1, sizeof( gsl_function ) );
    c->h = 1.0;
    c->omega_m = 0.25;
    c->omega_l = 1.0 - c->omega_m;
    c->w0 = -1.0;
    c->wa = 0.0;
}

void
cosmo_clean( COSMO * c )
{
    CHECK_FREE( c->df );
    gsl_integration_workspace_free( c->work );
}

COSMO *
cosmo_init( void )
{
    COSMO *c;
    c = ( COSMO * ) check_alloc( 1, sizeof( COSMO ) );
    cosmo_alloc( c );
    return c;
}

COSMO *
cosmo_kill( COSMO * c )
{
    cosmo_clean( c );
    CHECK_FREE( c );
    return NULL;
}

static inline double
pow2( const double x )
{
    return x * x;
}

static inline double
pow3( const double x )
{
    return x * x * x;
}

static inline double
pow4( const double x )
{
    return x * x * x * x;
}

static inline void
cosmo_set_omega_m( COSMO * c, const double om )
{
    c->omega_m = om;
}

static inline void
cosmo_set_omega_l( COSMO * c, const double ol )
{
    c->omega_l = ol;
}

static inline void
cosmo_set_flat_om( COSMO * c, const double omega_m )
{
    c->omega_m = omega_m;
    c->omega_l = 1.0 - omega_m;
}

static inline void
cosmo_set_h( COSMO * c, const double h )
{
    c->h = h;
}
static inline void
cosmo_set_w0( COSMO * c, const double w0)
{
    c->w0 = w0;
}
static inline void
cosmo_set_wa( COSMO * c, const double wa)
{
    c->wa = wa;
}

static inline double
cosmo_omega_m( COSMO * c )
{
    return c->omega_m;
}

static inline double
cosmo_omega_l( COSMO * c )
{
    return c->omega_l;
}

static inline double
cosmo_h( COSMO * c )
{
    return c->h;
}
static inline double
cosmo_w0( COSMO * c)
{
    return c->w0;
}
static inline double
cosmo_wa( COSMO * c)
{
    return c->wa;
}

static inline double
cosmo_chi( COSMO const *const c, const double z )
{
    double om, ol, ok, zp1,w0,wa;
    om = c->omega_m;
    ol = c->omega_l;
    w0 = c->w0;
    wa = c->wa;
    /*printf("w0wa floats: %4.2f   %4.2f", w0,wa); */
    ok = 1.0 - om - ol;
    zp1 = z + 1.0;
    return sqrt( om * pow3( zp1 ) + ok * pow2( zp1 ) + ol*powf(zp1,3*(1+w0+wa))*exp(-3*wa*(1-1/zp1))); /* *powf(zp1, 3*(1+w0+wa*(1-1/zp1)) ) ); */
}

/* cosmological distance integrand : CAN NOT INLINE! */
double
cosmo_cdi( double z, void *params )
{
    COSMO *c;
    c = ( COSMO * ) params;

    return 1.0 / cosmo_chi( c, z );
}

static inline double
cosmo_dist_hubble( COSMO const *const c )
{
    return ( SPEED_OF_LIGHT / ( c->h * HUBBLE_UNITS ) );
}

static inline double
cosmo_dist_co_los( COSMO const *const c, const double z )
{
    double res, err;
    c->df->function = &cosmo_cdi;
    c->df->params = ( void * ) c;

    gsl_integration_qag( c->df, 0.0, z, 0.0, 9.99999e-9,
                         COSMO_INT_SIZE, GSL_INTEG_GAUSS41, c->work, &res, &err );

    return cosmo_dist_hubble( c ) * res;
}

#endif
