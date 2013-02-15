#pragma once
#ifndef CT_INCLUDED
#define CT_INCLUDED

#include <math.h>

#ifndef PI
#define PI (4.0 * atan(1.0))
#define RAD2DEG (180.0 / PI)
#endif

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

static inline void
ct_rotate_about_zaxis( const double theta, double *const v )
{
    double x = v[0];
    double y = v[1];
    double st = sin( theta );
    double ct = cos( theta );

    v[0] = ct * x - st * y;
    v[1] = st * x + ct * y;
}

static inline void
ct_rotate_about_yaxis( const double theta, double *const v )
{
    double x = v[0];
    double z = v[2];
    double st = sin( theta );
    double ct = cos( theta );

    v[0] = ct * x - st * z;
    v[2] = st * x + ct * z;
}

static inline void
ct_rotate_about_xaxis( const double theta, double *const v )
{
    double y = v[1];
    double z = v[2];
    double st = sin( theta );
    double ct = cos( theta );

    v[1] = ct * y - st * z;
    v[2] = st * y + ct * z;
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
ct_xyz_to_sky( double const *const v, double *az, double *el )
{
    *el = acos( v[2] );
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
    ct_xyz_to_sky( v, ra, dec );
    *ra *= RAD2DEG;
    *dec *= RAD2DEG;
    *dec = 90.0 - *dec;
}

#endif
