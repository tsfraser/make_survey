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
 *   1. clean c-only code without globals.
 *   2. few to no dependencies.
 *   3. simple to use and utilize within larger codes.
 *
 *   There is no attempt to provide the complete functionality
 *   that MANGLE does, nor have flexible input.  Please use the
 *   full MANGLE utilties if you need to convert to an accepted
 *   polygon format that this code can handle.
 *
 *   Prefix for functions is "mply_", short for mangle polygon.
 *
 *   NOTE: there is no attempt to separate "internal" and "external"
 *   functions by name.
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
#ifdef M_PI
#define PI M_PI
#else
#define PI (4.0 * atan(1.0))
#endif
#endif

#define DEG2RAD ( PI / 180.0 )

typedef int MANGLE_INT;         /* signed integer */

typedef struct {
    double x[3];
} MANGLE_VEC;

typedef struct {
    double x[3];
    double m;
} MANGLE_CAP;

typedef struct {
    MANGLE_INT ipoly;           /* internal index: will be unique */
    MANGLE_INT polyid;
    MANGLE_INT pixel;
    MANGLE_INT ncap;
    MANGLE_CAP *cap;
    double weight;
    double area;
} MANGLE_POLY;

/* this is a linked list structure */
typedef struct {
    void *data;
    void *next;
} DATA_LIST;

typedef struct {
    MANGLE_INT npoly;
    MANGLE_POLY *poly;
    MANGLE_INT pix_res;         /* pix_res = 0 is full sky: aka no pixels */
    DATA_LIST *pix;             /* pixel-indexed array of linked-lists */
} MANGLE_PLY;

/* this is a utility function for simple pixelization scheme */
static inline size_t
mply_pow2i( const int x )
{
    /* as long as x is small, do pow() via bitshift! */
    return ( 1 << x );
}

static inline size_t
mply_pix_count( const int res )
{
    size_t npix;
    if( res < 1 )
        return 0;
    npix = mply_pow2i( 2 * res );
    return npix;
}

void
mply_pix_alloc( MANGLE_PLY * const ply, const int pix_res )
{
    /* use pix_res to allocate ply->pix array */
    size_t i, count;
    count = mply_pix_count( pix_res );
    ply->pix = ( DATA_LIST * ) check_alloc( count, sizeof( DATA_LIST ) );
    ply->pix_res = pix_res;

    /* for sanity, ensure all starting points are NULL */
    for( i = 0; i < count; i++ ) {
        ply->pix[i].data = NULL;        /* pointer to MANGLE_POLY */
        ply->pix[i].next = NULL;        /* pointer to next DATA_LIST */
    }
}

void
mply_pix_clean( MANGLE_PLY * const ply )
{
    size_t i, count;
    DATA_LIST *cl;              /* current list */
    DATA_LIST *nl;              /* next list */

    count = mply_pix_count( ply->pix_res );

    /* pixels stored as a linked-list, have to decend to free() */
    for( i = 0; i < count; i++ ) {
        cl = &( ply->pix[i] );
        nl = ( DATA_LIST * ) cl->next;
        /* NOTE: this skips the first element, which is INTENDED! */
        while( nl != NULL ) {
            cl = nl;
            nl = ( DATA_LIST * ) cl->next;
            CHECK_FREE( cl );
        }
    }
    /* now we can clean up the original array */
    CHECK_FREE( ply->pix );
    ply->pix_res = 0;
}

/* next several routines modeled after code in which_pixel.c in original mangle code */
/* Given pixel resolution, what is the ID of the starting pixel? */
static inline MANGLE_INT
mply_pix_id_start( MANGLE_PLY const *const ply )
{
    int pix_id;
    MANGLE_INT res;

    res = ply->pix_res;
    pix_id = ( mply_pow2i( 2 * res ) - 1 ) / 3;
    return pix_id;
}

/* INDEX refers to the internal storage index, which is in pixel order but
 * zero-indexed rather than numbered according to resolution as the
 * "simple pixelization" scheme in MANGLE does */
static inline MANGLE_INT
mply_pix_which_index( MANGLE_PLY const *const ply, const double az, double el )
{
    int n, m;
    MANGLE_INT base_pix, pow2r;

    pow2r = mply_pow2i( ply->pix_res );

    /* algorithm made to replicate comparisons in Mangle's which_pixel.c */
    if( sin( el ) == 1.0 ) {
        n = 0;
    } else {
        n = ( int ) ceil( ( 1.0 - sin( el ) ) / 2.0 * pow2r ) - 1;
    }
    m = ( int ) floor( az / 2.0 / PI * pow2r );
    base_pix = pow2r * n + m;

    return base_pix;
}

static inline MANGLE_INT
mply_pix_which_id( MANGLE_PLY const *const ply, const double az, const double el )
{
    return mply_pix_which_index( ply, az, el ) + mply_pix_id_start( ply );
}

static inline MANGLE_INT
mply_pix_index_from_id( MANGLE_PLY const *const ply, MANGLE_INT id )
{
    MANGLE_INT index;
    size_t npix = mply_pix_count( ply->pix_res );
    index = id - mply_pix_id_start( ply );

    if( index < 0 || ( size_t ) index >= npix ) {
        fprintf( stderr,
                 "MANGLE Error: pixel_id=%zd does not match with pixel_res=%zd.\n",
                 ( ssize_t ) id, ( ssize_t ) ply->pix_res );
        exit( EXIT_FAILURE );
    }
    return index;
}

static inline size_t
mply_pix_npoly( MANGLE_PLY const *const ply, MANGLE_INT ipix )
{
    size_t count = 0;
    DATA_LIST *dl;

    dl = ( DATA_LIST * ) & ply->pix[ipix];
    while( dl != NULL ) {
        if( dl->data != NULL ) {
            count += 1;
        }
        dl = ( DATA_LIST * ) dl->next;
    }
    return count;
}

void
mply_pix_addpoly( MANGLE_PLY const *const ply, MANGLE_POLY const *const p )
{
    MANGLE_INT index;
    DATA_LIST *dl;
    if( ply->pix_res < 1 ) {
        fprintf( stderr,
                 "MANGLE Error: Tried to add pixel without proper PIXEL initialization!\n" );
        exit( EXIT_FAILURE );
    }
    index = mply_pix_index_from_id( ply, p->pixel );
    dl = ( DATA_LIST * ) & ply->pix[index];
    while( dl->data != NULL ) {
        if( dl->next != NULL ) {
            dl = ( DATA_LIST * ) dl->next;
        } else {
            dl->next = ( void * ) check_alloc( 1, sizeof( DATA_LIST ) );
            dl = ( DATA_LIST * ) dl->next;
            /* values of this new structure get set below */
            break;
        }
    }
    dl->data = ( void * ) p;    /* pointer to polygon structure */
    dl->next = NULL;            /* pointer to next DATA_LIST */
}

/* The POLY structure holds a list of caps, not the umbrella PLY structure */
void
mply_poly_alloc( MANGLE_POLY * p, const MANGLE_INT ipoly, const MANGLE_INT polyid,
                 const MANGLE_INT ncap, const double weight, const MANGLE_INT pixel,
                 const double area )
{
    p->ipoly = ipoly;
    p->polyid = polyid;
    p->cap = ( MANGLE_CAP * ) check_alloc( ncap, sizeof( MANGLE_CAP ) );
    p->ncap = ncap;
    p->weight = weight;
    p->pixel = pixel;
    p->area = area;
}

void
mply_poly_clean( MANGLE_POLY * p )
{
    p->polyid = -1;
    CHECK_FREE( p->cap );
    p->ncap = 0;
    p->weight = 0.0;
    p->pixel = -1;
    p->area = 0.0;
}

/* This is the main polygon structure */
void
mply_alloc( MANGLE_PLY * ply, MANGLE_INT npoly )
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
mply_clean( MANGLE_PLY * ply )
{
    MANGLE_INT i;
    MANGLE_POLY *p;
    for( i = 0; i < ply->npoly; i++ ) {
        p = &( ply->poly[i] );
        mply_poly_clean( p );
    }
    CHECK_FREE( ply->poly );
    ply->npoly = 0;
    if( ply->pix_res > 0 )
        mply_pix_clean( ply );
}

MANGLE_PLY *
mply_init( MANGLE_INT npoly )
{
    MANGLE_PLY *ply;
    ply = ( MANGLE_PLY * ) check_alloc( 1, sizeof( MANGLE_PLY ) );
    mply_alloc( ply, npoly );
    return ply;
}

MANGLE_PLY *
mply_kill( MANGLE_PLY * ply )
{
    mply_clean( ply );
    CHECK_FREE( ply );
    return NULL;
}

void
mply_read_file_into( MANGLE_PLY * const ply, char const *const filename )
{
    /* read in polygon format */
    int check;
    int npoly = 0;
    MANGLE_INT ipoly = 0;
    simple_reader *sr;
    char *line;

    mply_clean( ply );

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

    mply_alloc( ply, npoly );

    /* read other header directives */
    while( sr_readline( sr ) ) {
        if( sr_line_isempty( sr ) )
            continue;

        line = sr_line( sr );

        if( strncmp( "polygon", line, 7 ) == 0 ) {
            /* we're done reading headers */
            break;
        }

        if( strncmp( "pixelization", line, 12 ) == 0 ) {
            int res;
            check = sscanf( line, "pixelization %ds", &res );
            if( check != 1 ) {
                fprintf( stderr,
                         "MANGLE Warning: Only simple pixel scheme is currently supported: %s\n",
                         sr_filename( sr ) );
                continue;
            }
            mply_pix_alloc( ply, res );
        }
    }

    /* now off to POLY processing */
    ipoly = 0;
    do {
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
            mply_poly_alloc( p, ipoly, polyid, ncap, weight, pixel, area );
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
            if( ply->pix_res > 0 ) {
                mply_pix_addpoly( ply, p );
            }
            ipoly += 1;
        }
        /* silently ignore anything else, only processing polygons here */
    } while( sr_readline( sr ) );

    sr_kill( sr );

    if( ipoly != ply->npoly ) {
        fprintf( stderr,
                 "MANGLE Error: bad number of polygons read! Expected %zd, read %zd\n",
                 ( ssize_t ) ply->npoly, ( ssize_t ) ipoly );
        exit( EXIT_FAILURE );
    }
}

MANGLE_PLY *
mply_read_file( char const *const filename )
{
    MANGLE_PLY *ply;
    ply = mply_init( 0 );
    mply_read_file_into( ply, filename );
    return ply;
}

/* this can be abstracted: a calling code can just use (void *) */
MANGLE_VEC *
mply_vec_init( void )
{
    MANGLE_VEC *vec3;
    vec3 = ( MANGLE_VEC * ) check_alloc( 1, sizeof( MANGLE_VEC ) );
    return vec3;
}

MANGLE_VEC *
mply_vec_kill( MANGLE_VEC * vec3 )
{
    CHECK_FREE( vec3 );
    return vec3;
}

/* convert polor sky coordinates to unit vector:
 * el = elevation / polar
 * az = azimuthal
 */
static inline void
mply_vec_from_polar( MANGLE_VEC * vec3, const double az, const double el )
{
    /* simplified mapping to match MANGLE internals */
    vec3->x[0] = cos( el ) * cos( az );
    vec3->x[1] = cos( el ) * sin( az );
    vec3->x[2] = sin( el );
}

static inline void
mply_vec_from_radec( MANGLE_VEC * vec3, const double ra, const double dec )
{
    /* simplified mapping to match MANGLE internals */
    mply_vec_from_polar( vec3, ra * DEG2RAD, dec * DEG2RAD );
}

static inline MANGLE_INT
mply_within_cap( MANGLE_CAP const *const cap, MANGLE_VEC const *const vec3 )
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
mply_within_poly( MANGLE_POLY const *const p, MANGLE_VEC const *const vec3 )
{
    MANGLE_INT i;
    MANGLE_CAP *c;
    c = p->cap;
    for( i = 0; i < p->ncap; i++ ) {
        if( !mply_within_cap( &c[i], vec3 ) )
            return FALSE;
    }
    return TRUE;
}

/* short circuit: finds FIRST matching polygon and does not continue checking! */
static inline MANGLE_INT
mply_find_polyindex_vec( MANGLE_PLY const *const ply, MANGLE_VEC const *const vec3 )
{
    MANGLE_INT i;
    MANGLE_POLY *p;

    for( i = 0; i < ply->npoly; i++ ) {
        p = &( ply->poly[i] );
        if( mply_within_poly( p, vec3 ) )
            return i;
    }
    return -1;
}

/* In retrospect, the pixel code is probably unnecessarily complex.
 *
 * But hey, it's working now... so why bother changing it.
 *
 * What would have probably worked fine, and avoided the whole linked-list thing,
 * is just to loop over all polygons (brute force) and not check any caps unless
 * the pixel_id corresponded to the ra/dec pair. Blindingly obvious now, but it
 * didn't occur to me right away.
 *
 * The current method avoids having to do such a large loop (more efficient in
 * that sense). However, it has several drawbacks:
 * 1. it requires more memory (pix list and all linked-lists elements
 * 2. iterates over linked-list, so more random access memory
 */
static inline MANGLE_INT
mply_find_polyindex_pix( MANGLE_PLY const *const ply, const double az, const double el )
{
    MANGLE_INT ipix;
    MANGLE_POLY *p;
    MANGLE_VEC vec3;
    DATA_LIST *dl;

    mply_vec_from_polar( &vec3, az, el );
    ipix = mply_pix_which_index( ply, az, el );
    dl = ( DATA_LIST * ) & ply->pix[ipix];

    /* transverse our linked-list and test for matches */
    while( dl != NULL && dl->data != NULL ) {
        p = ( MANGLE_POLY * ) dl->data;
        dl = ( DATA_LIST * ) dl->next;
        if( mply_within_poly( p, &vec3 ) )
            return p->ipoly;
    }

    return -1;
}

static inline MANGLE_INT
mply_find_polyindex_polar( MANGLE_PLY const *const ply, const double az, const double el )
{
    MANGLE_INT i = -1;

    if( ply->pix_res > 0 ) {
        i = mply_find_polyindex_pix( ply, az, el );
    } else {
        MANGLE_VEC vec3;
        mply_vec_from_polar( &vec3, az, el );
        i = mply_find_polyindex_vec( ply, &vec3 );
    }

    return i;
}

static inline MANGLE_INT
mply_find_polyindex_radec( MANGLE_PLY const *const ply, const double ra, const double dec )
{
    return mply_find_polyindex_polar( ply, ra * DEG2RAD, dec * DEG2RAD );
}

MANGLE_POLY *
mply_poly_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    if( index >= ply->npoly || index < 0 ) {
        fprintf( stderr, "MANGLE Error: invalid POLY index: %zd\n", ( ssize_t ) index );
        exit( EXIT_FAILURE );
    }

    return &( ply->poly[index] );
}

static inline MANGLE_INT
mply_polyid_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    MANGLE_POLY *p;
    if( index < 0 )
        return -1;
    p = mply_poly_from_index( ply, index );

    return p->polyid;
}

static inline double
mply_weight_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    MANGLE_POLY *p;
    if( index < 0 )
        return 0.0;
    p = mply_poly_from_index( ply, index );

    return p->weight;
}

static inline double
mply_area_from_index( MANGLE_PLY const *const ply, const MANGLE_INT index )
{
    MANGLE_POLY *p;
    if( index < 0 )
        return 0.0;
    p = mply_poly_from_index( ply, index );

    return p->area;
}

static inline double
mply_area_total( MANGLE_PLY const *const ply, const double min_weight )
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
mply_area_weighted_total( MANGLE_PLY const *const ply, const double min_weight )
{
    MANGLE_INT i;
    double wa = 0.0;
    for( i = 0; i < ply->npoly; i++ ) {
        if( ply->poly[i].weight < min_weight )
            continue;
        wa += ( ply->poly[i].area * ply->poly[i].weight );
    }
    return wa;
}

static inline MANGLE_INT
mply_find_polyid( MANGLE_PLY const *const ply, MANGLE_VEC const *const vec3 )
{
    MANGLE_INT index = mply_find_polyindex_vec( ply, vec3 );

    if( index < 0 )
        return -1;

    return mply_polyid_from_index( ply, index );
}

#endif
