#pragma once
#ifndef CT_INCLUDED
#define CT_INCLUDED

#include <math.h>

#ifndef PI
#ifdef M_PI
#define PI M_PI
#else
#define PI (4.0 * atan(1.0))
#endif
#endif

#define RAD2DEG (180.0 / PI)

static inline double
ct_rad_from_deg( double const deg )
{
    return deg / RAD2DEG;
}

static inline double
ct_deg_from_rad( double const rad )
{
    return rad * RAD2DEG;
}

/* below should generalize the other axis rotations */
static inline void
ct_rotate_about_axis( const double theta, double *const v, const int axis )
{
    int i, ia, ic[2];
    double a, b;
    double st = sin( theta );
    double ct = cos( theta );

    /* select which axes to manipulate:
     * basically, not the requested one in the order of X,Y,Z */
    ia = 0;
    for( i = 0; i < 3; i++ ) {
        if( i != axis ) {
            ic[ia] = i;
            ia += 1;
        }
    }

    a = v[ic[0]];
    b = v[ic[1]];

    v[ic[0]] = ct * a - st * b;
    v[ic[1]] = st * a + ct * b;
}

static inline void
ct_rotate_about_z( const double theta, double *const v )
{
    ct_rotate_about_axis( theta, v, 2 );
}

static inline void
ct_rotate_about_y( const double theta, double *const v )
{
    ct_rotate_about_axis( theta, v, 1 );
}

static inline void
ct_rotate_about_x( const double theta, double *const v )
{
    ct_rotate_about_axis( theta, v, 0 );
}

static inline double
ct_vec_dot( double const *const v1, double const *const v2 )
{
    int i;
    double r = 0.0;

    for( i = 0; i < 3; i++ )
        r += v1[i] * v2[i];

    return r;
}

static inline double
ct_vec_mag( double const *const v )
{
    return sqrt( ct_vec_dot( v, v ) );
}

static inline void
ct_xyz_to_polar( double const *const v, double *az, double *el )
{
    double r = ct_vec_mag( v );
    *el = acos( v[2] / r );
    /* algorithm below for compatablility with some other codes */
    *az = atan( v[1] / v[0] );
    if( v[0] < 0 )
        *az += PI;
    else if( v[1] < 0 )
        *az += 2.0 * PI;
    else if( v[0] == 0 && v[1] == 0 )
        *az = 0.0;
}

static inline void
ct_xyz_to_radec( double const *const v, double *ra, double *dec )
{
    ct_xyz_to_polar( v, ra, dec );
    *ra *= RAD2DEG;
    *dec *= RAD2DEG;
    *dec = 90.0 - *dec;
}

static inline void
ct_polar_to_xyz( double const az, double const el, double *const x )
{
    x[0] = sin( el ) * cos( az );
    x[1] = sin( el ) * sin( az );
    x[2] = cos( el );
}

static inline void
ct_radec_to_xyz( double const ra, double const dec, double *const x )
{
    double az = ra / RAD2DEG;
    double el = ( 90 - dec ) / RAD2DEG;
    ct_polar_to_xyz( az, el, x );
}

#endif
