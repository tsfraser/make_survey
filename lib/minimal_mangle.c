/* simplified C-only code utilities to minimally use mangle polygons
 *
 * Cameron McBride
 * cameron.mcbride AT gmail.com
 * January 2013
 */
#pragma once
#ifndef MINIMAL_MANGLE_INCLUDED
#define MINIMAL_MANGLE_INCLUDED

/*
 * Some description:
 *   The main focus of these utilites:
 *   1. clean c-only code
 *   2. few dependencies
 *   3. simple use and easy to utilize within larger codes
 *
 *   There is no attempt to provide the complete functionality
 *   that MANGLE does, nor have near as flexible input.  Please use
 *   the full MANGLE utilties if you need to convert to an accepted
 *   polygon format that this code can read.
 **/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <check_alloc.c>
#include <check_fopen.c>
#include <simple_reader.c>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#ifndef PI
#define PI (4.0 * atan(1.0))
#endif

typedef int MANGLE_INT;         /* signed integer */

typedef struct {
    double x[3];
} MANGLE_VEC;

typedef struct {
    double x[3];
    double m;
} MANGLE_CAP;

typedef struct {
    MANGLE_INT polyid;
    MANGLE_INT pixel;
    MANGLE_INT ncap;
    MANGLE_CAP *cap;
    double weight;
    double area;
} MANGLE_POLY;

typedef struct {
    void *data;
    void *next;
} DATA_LIST;

typedef struct {
    MANGLE_INT npoly;
    MANGLE_POLY *poly;
    MANGLE_INT pix_res;
    DATA_LIST *pix;             /* pixel-indexed array of linked-lists */
} MANGLE_PLY;

enum {
    PLY_PIX_NONE = 0,
    PLY_PIX_SIMPLE = 1,
    PLY_PIX_SDSSPIX = 2
};

void
ply_poly_alloc( MANGLE_POLY * p, const MANGLE_INT polyid, const MANGLE_INT ncap,
                const double weight, const MANGLE_INT pixel, const double area )
{
    p->polyid = polyid;
    p->cap = ( MANGLE_CAP * ) check_alloc( ncap, sizeof( MANGLE_CAP ) );
    p->ncap = ncap;
    p->weight = weight;
    p->pixel = pixel;
    p->area = area;
}

void
ply_poly_clean( MANGLE_POLY * p )
{
    p->polyid = -1;
    CHECK_FREE( p->cap );
    p->ncap = 0;
    p->weight = 0.0;
    p->pixel = -1;
    p->area = 0.0;
}

void
ply_alloc( MANGLE_PLY * ply, MANGLE_INT npoly )
{
    if( npoly > 0 ) {
        ply->poly = ( MANGLE_POLY * ) check_alloc( npoly, sizeof( MANGLE_POLY ) );
    } else {
        ply->poly = NULL;
    }
    ply->npoly = npoly;
    ply->pix_res = 0;
}

void
ply_clean( MANGLE_PLY * ply )
{
    MANGLE_INT i;
    MANGLE_POLY *p;
    for( i = 0; i < ply->npoly; i++ ) {
        p = &( ply->poly[i] );
        ply_poly_clean( p );
    }
    CHECK_FREE( ply->poly );
    ply->npoly = 0;
    ply->pix_res = 0;
}

MANGLE_PLY *
ply_init( MANGLE_INT npoly )
{
    MANGLE_PLY *ply;
    ply = ( MANGLE_PLY * ) check_alloc( 1, sizeof( MANGLE_PLY ) );
    ply_alloc( ply, npoly );
    return ply;
}

MANGLE_PLY *
ply_kill( MANGLE_PLY * ply )
{
    ply_clean( ply );
    CHECK_FREE( ply );
    return NULL;
}

void
ply_read_file_into( MANGLE_PLY * const ply, char const *const filename )
{
    /* read in polygon format */
    int check;
    int npoly = 0;
    MANGLE_INT ipoly = 0;
    simple_reader *sr;
    char *line;

    ply_clean( ply );

    /* read-by-line and process MANGLE_PLY format file */
    sr = sr_init( filename );

    /* first line sets up the polygons */
    line = sr_readline( sr );
    check = sscanf( line, "%d polygons", &npoly );
    if( check != 1 || npoly < 1 ) {
        fprintf( stderr,
                 "MANGLE Error: polygons (%d) must be positive in file: %s\n",
                 npoly, sr_filename( sr ) );
        exit( EXIT_FAILURE );
    }

    ply_alloc( ply, npoly );

    ipoly = 0;
    while( sr_readline( sr ) ) {
        int i, polyid, ncap, pixel;
        double weight, area;
        MANGLE_POLY *p;

        if( sr_line_isempty( sr ) )
            continue;

        line = sr_line( sr );

        if( strncmp( "polygon", line, 7 ) == 0 ) {
            check =
                sscanf( line,
                        "polygon %d ( %d caps, %lf weight, %d pixel, %lf",
                        &polyid, &ncap, &weight, &pixel, &area );
            if( check != 5 || ncap < 1 ) {
                fprintf( stderr,
                         "MANGLE Error: polygon read error line %d in file: %s\n",
                         sr_linenum( sr ), sr_filename( sr ) );
                exit( EXIT_FAILURE );
            }

            if( ipoly >= ply->npoly ) {
                fprintf( stderr,
                         "MANGLE Error: too many polygons on line %d in file: %s\n",
                         sr_linenum( sr ), sr_filename( sr ) );
                exit( EXIT_FAILURE );
            }

            /* we're starting a valid polygon! */
            p = &ply->poly[ipoly];
            ply_poly_alloc( p, polyid, ncap, weight, pixel, area );
            for( i = 0; i < ncap; i++ ) {
                MANGLE_CAP *c;
                c = &p->cap[i];
                line = sr_readline( sr );
                check = sscanf( line, "%lf %lf %lf %lf", &c->x[0], &c->x[1], &c->x[2], &c->m );
                if( check != 4 ) {
                    fprintf( stderr,
                             "MANGLE Error: cap read error on line %d in file: %s\n",
                             sr_linenum( sr ), sr_filename( sr ) );
                    exit( EXIT_FAILURE );
                }
            }
            ipoly += 1;
        }
        /* XXX pretty much ignore anything else, but will implement simple pixelization */
    }

    sr_kill( sr );

    if( ipoly != ply->npoly ) {
        fprintf( stderr,
                 "MANGLE Error: too few polygons read! Expected %zd, read %zd\n",
                 ( ssize_t ) ply->npoly, ( ssize_t ) ipoly );
        exit( EXIT_FAILURE );
    }
}

MANGLE_PLY *
ply_read_file( char const *const filename )
{
    MANGLE_PLY *ply;
    ply = ply_init( 0 );
    ply_read_file_into( ply, filename );
    return ply;
}

/* this can be abstracted: a calling code can just use (void *) */
MANGLE_VEC *
ply_vec_init( void )
{
    MANGLE_VEC *vec3;
    vec3 = ( MANGLE_VEC * ) check_alloc( 1, sizeof( MANGLE_VEC ) );
    return vec3;
}

MANGLE_VEC *
ply_vec_kill( MANGLE_VEC * vec3 )
{
    CHECK_FREE( vec3 );
    return vec3;
}

/* convert polor sky coordinates to unit vector:
 * el = elevation / polar
 * az = azimuthal
 */
static inline void
ply_vec_from_polar( MANGLE_VEC * vec3, const double az, const double el )
{
    vec3->x[0] = sin( el ) * cos( az );
    vec3->x[1] = sin( el ) * sin( az );
    vec3->x[2] = cos( el );
}

static inline void
ply_vec_from_radec_deg( MANGLE_VEC * vec3, const double ra, const double dec )
{
    double az = PI / 180.0 * ra;
    double el = PI / 180.0 * ( 90.0 - dec );
    ply_vec_from_polar( vec3, az, el );
}

static inline MANGLE_INT
ply_within_cap( MANGLE_CAP const *const cap, MANGLE_VEC const *const vec3 )
{
    const double *c;
    const double *v;
    c = cap->x;
    v = vec3->x;
    double cd = 1.0 - c[0] * v[0] - c[1] * v[1] - c[2] * v[2];
    if( cap->m < 0.0 ) {
        if( cd > fabs( cap->m ) )
            return TRUE;
    } else {
        if( cd < cap->m )
            return TRUE;
    }

    return FALSE;
}

static inline MANGLE_INT
ply_within_poly( MANGLE_POLY const *const p, MANGLE_VEC const *const vec3 )
{
    MANGLE_INT i;
    MANGLE_CAP *c;
    c = p->cap;
    for( i = 0; i < p->ncap; i++ ) {
        if( !ply_within_cap( &c[i], vec3 ) )
            return FALSE;
    }
    return TRUE;
}

/* short circuit, finds *FIRST* matching polygon and does not continue checking! */
static inline MANGLE_INT
ply_find_index( MANGLE_PLY const *const ply, MANGLE_VEC const *const vec3 )
{
    MANGLE_INT i;
    MANGLE_POLY *p;

    for( i = 0; i < ply->npoly; i++ ) {
        p = &( ply->poly[i] );
        if( ply_within_poly( p, vec3 ) )
            return i;
    }
    return -1;
}

MANGLE_POLY *
ply_poly_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    if( index >= ply->npoly || index < 0 ) {
        fprintf( stderr, "MANGLE Error: invalid POLY index: %zd\n", ( ssize_t ) index );
        exit( EXIT_FAILURE );
    }

    return &( ply->poly[index] );
}

static inline MANGLE_INT
ply_polyid_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    return ( ply_poly_from_index( ply, index ) )->polyid;
}

static inline double
ply_weight_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    MANGLE_POLY *p;
    p = ply_poly_from_index( ply, index );
    return p->weight;
}

static inline double
ply_area_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    return ( ply_poly_from_index( ply, index ) )->area;
}

static inline double
ply_area_total( MANGLE_PLY const *const ply, const double min_weight )
{
    MANGLE_INT i;
    double area = 0.0;
    for( i = 0; i < ply->npoly; i++ ) {
        if( ply->poly[i].weight < min_weight )
            continue;
        area += ply->poly[i].area;
    }
    return area;
}

static inline double
ply_area_weighted_total( MANGLE_PLY const *const ply, const double min_weight )
{
    MANGLE_INT i;
    double warea = 0.0;
    for( i = 0; i < ply->npoly; i++ ) {
        if( ply->poly[i].weight < min_weight )
            continue;
        warea += ( ply->poly[i].area * ply->poly[i].weight );
    }
    return warea;
}

static inline MANGLE_INT
ply_find_polyid( MANGLE_PLY const *const ply, MANGLE_VEC const *const vec3 )
{
    MANGLE_INT index = ply_find_index( ply, vec3 );

    if( index < 0 )
        return -1;

    return ply_polyid_from_index( ply, index );
}

/* PIXEL routines: NOT YET TESTED .. STILL IN DEVELOPMENT.  LIKELY BROKEN!! */
static inline int
ply_pow2i( const int x )
{
    /* as long as x is small, do pow() via bitshift! */
    return ( 1 << x );
}

static inline size_t
ply_pix_count( const int res )
{
    size_t npix;
    if( res < 1 )
        return 0;
    npix = ply_pow2i( 2 * res );
    return npix;
}

void
ply_pix_init( MANGLE_PLY * const ply, const int pix_res )
{
    /* use pix_res to allocate ply->pix array */
    size_t count = ply_pix_count( pix_res );
    ply->pix = ( DATA_LIST * ) check_alloc( count, sizeof( DATA_LIST ) );
    ply->pix_res = pix_res;
}

void
ply_pix_addpoly( MANGLE_PLY const *const ply, MANGLE_POLY const *const p )
{
    // XXX in development
    // XXX check PIXEL_ID matches the implied resolution!
    if( ply->pix_res < 1 ) {
        fprintf( stderr,
                 "Error: Tried to add pixel without proper initializing PIXEL structure!\n" );
    }
}

/* next several routines modeled after code in which_pixel.c in original mangle code */

/* Given pixel resolution, what is the ID of the starting pixel? */
static inline MANGLE_INT
ply_pix_id_start( MANGLE_PLY const *const ply )
{
    int pix_id;
    MANGLE_INT res;

    res = ply->pix_res;
    pix_id = ( ply_pow2i( 2 * res ) - 1 ) / 3;
    return pix_id;
}

/* INDEX refers to the internal storage index, which is in pixel order but
 * zero-indexed rather than numbered according to resolution as the
 * "simple pixelization" scheme in MANGLE does */
static inline MANGLE_INT
ply_pix_which_index( MANGLE_PLY const *const ply, const double az, const double el )
{
    int n, m, base_pix;
    MANGLE_INT pow2r;

    pow2r = ply_pow2i( ply->pix_res );

    if( fabs( sin( el ) - 1.0 ) < 1e-10 ) {
        n = 0;
    } else {
        n = ceil( ( 1.0 - sin( el ) ) / 2.0 * pow2r ) - 1;
    }
    m = floor( az / 2.0 / PI );
    base_pix = ( n + m ) * pow2r;

    return base_pix;
}

static inline MANGLE_INT
ply_pix_which_id( MANGLE_PLY const *const ply, const double az, const double el )
{
    return ply_pix_which_index( ply, az, el ) + ply_pix_id_start( ply );
}

#endif
