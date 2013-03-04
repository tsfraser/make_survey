/* make_survey
 * -----------
 *
 * GOAL: take points distributed in a periodic box (mock galaxies, halos, etc)
 * and optionally apply various transformations to make a mock survey. The basic steps
 * that can be applied are:
 * 1. BoxRemap
 * 2. Translate / Rotate
 * 3. Apply redshift distortions
 * 4. Sky projection (using proper cosmological distance)
 * 5. Trim to survey footprint (mangle mask: ascii PLY)
 * 6. Downsample based on sky completeness (in polygon mask)
 * 7. Downsample based on redshift completeness (additional input)
 *
 *
 * This is glue code, with most functionality in the "library" codes.  Most is written
 * in C, but this is a modified version of BoxRemap C++ version from Carlson & White:
 * http://mwhite.berkeley.edu/BoxRemap/
 *
 *
 * Cameron McBride
 * cameron.mcbride@gmail.com
 * February 2013
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <cuboid.h>

#include <simple_reader.c>
#include <check_fopen.c>

#include <config.c>
#include <rng_gsl.c>
#include <spline_gsl.c>
#include <coord_transforms.c>
#include <cosmology.c>
#include <zsel.c>
#include <minimal_mangle.c>

#define NSPLINE 1000            /* number of elements to train cosmology spline with */

#define ON 1
#define OFF 0

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

int
main( int argc, char *argv[] )
{
    char *file_config = NULL;
    char *file_mock = NULL;
    char *file_out = NULL;
    FILE *fdout;

    size_t nread = 0, nout = 0;

    char *line;
    simple_reader *sr;

    CONFIG *conf;
    void *rng;
    ZSEL *zsel;
    COSMO *cosmo;
    SPLINE *spl;
    MANGLE_PLY *ply;
    double rmin, rmax;

    ply = NULL;

    if( argc != 4 ) {
        fprintf( stderr, "USAGE:\t%s  CONFIG_FILE  MOCK_IN  RDZ_OUT\n", argv[0] );
        exit( EXIT_FAILURE );
    }

    file_config = argv[1];
    file_mock = argv[2];
    file_out = argv[3];         /* XXX generate this!! */
    if( strcmp( file_mock, file_out ) == 0 ) {
        fprintf( stderr, "Error: OUTPUT file must be different than INPUT!\n" );
        exit( EXIT_FAILURE );
    }
    // XXX read config file!
    conf = conf_init(  );
    {
        conf->file_mask = file_config;  // XXX hack!!
        conf->downsample_sky = ON;

        conf->lbox = 2560.0;
        conf->omega_m = 0.29;
        conf->omega_l = 1.0 - conf->omega_m;
        conf->hubble = 1.0;
        conf->seed = 6452;

        conf->min_sky_weight = 0.7;
        conf->zspace = 1;
        conf->zin = 0.55;

        conf->zmin = 0.40;
        conf->zmax = 0.75;

        conf->u[0] = 1;
        conf->u[1] = 1;
        conf->u[2] = 0;

        conf->u[3] = 0;
        conf->u[4] = 0;
        conf->u[5] = 1;

        conf->u[6] = 1;
        conf->u[7] = 0;
        conf->u[8] = 0;

        /* NGC */
        conf->pre_rot[0] = 1;
        conf->pre_rot[1] = 0;
        conf->pre_rot[2] = 0;

        conf->t[0] = -1810.0;
        conf->t[1] = -350.0;
        conf->t[2] = -10.0;

        conf->rot[0] = 355.0;
        conf->rot[1] = 0.0;
        conf->rot[2] = 95.0;
        conf->file_zsel = "cmass-dr10v5-ngc.zsel";

        /* SGC */
//         conf->pre_rot[0] = 1;
//         conf->pre_rot[1] = 0;
//         conf->pre_rot[2] = 0;

//         conf->t[0] = -1810.0;
//         conf->t[1] = -310.0;
//         conf->t[2] = -380.0;

//         conf->rot[0] = 340.0;
//         conf->rot[1] = 0.0;
//         conf->rot[2] = 270.0;
//         conf->file_zsel = "cmass-dr10v5-sgc.zsel";

    }

    {
        fprintf( stderr, "pre-rotation> about each axis in 90 deg units [%d,%d,%d]\n",
                 conf->pre_rot[0], conf->pre_rot[1], conf->pre_rot[2] );
    }

    /* initialize cuboid remapping */
    Cuboid R( conf->u );
    {
        double r[3];
        // print out lattice vectors
        fprintf( stderr, "remap> LATTICE VECTORS:\n" );
        for( int i = 0; i < 3; i++ ) {
            fprintf( stderr, "remap>   u%d = (", i + 1 );
            for( int j = 0; j < 3; j++ ) {
                fprintf( stderr, " %d ", conf->u[i * 3 + j] );
                if( j < 2 )
                    fprintf( stderr, "," );
            }
            fprintf( stderr, ")\n" );
        }

        r[0] = R.L1 * conf->lbox;
        r[1] = R.L2 * conf->lbox;
        r[2] = R.L3 * conf->lbox;

        fprintf( stderr, "remap> REMAPPING max x: %8.2f to %6.2f\n", conf->lbox, r[0] );
        fprintf( stderr, "remap> REMAPPING max y: %8.2f to %6.2f\n", conf->lbox, r[1] );
        fprintf( stderr, "remap> REMAPPING max z: %8.2f to %6.2f\n", conf->lbox, r[2] );

        fprintf( stderr, "trans> TRANSLATE x: %8.2f to %6.2f\n", conf->t[0], r[0] + conf->t[0] );
        fprintf( stderr, "trans> TRANSLATE y: %8.2f to %6.2f\n", conf->t[1], r[1] + conf->t[1] );
        fprintf( stderr, "trans> TRANSLATE z: %8.2f to %6.2f\n", conf->t[2], r[2] + conf->t[2] );

        fprintf( stderr, "rot> ROTATION about x,y,z in degrees: %g,%g,%g\n",
                 conf->rot[0], conf->rot[1], conf->rot[2] );
    }

    /* do all initialization */
    {
        fprintf( stderr, "rng> initializing random number generator with seed = %u\n", conf->seed );
        rng = rng_init( conf->seed );

        /* initialize cosmology */
        cosmo = cosmo_init(  );
        cosmo_set_omega_m( cosmo, conf->omega_m );
        cosmo_set_omega_l( cosmo, conf->omega_l );
        cosmo_set_h( cosmo, 1.0 );

        fprintf( stderr, "cosmo> Cosmology: (Om, Ol, h) = (%g, %g, %.1f)\n",
                 cosmo_omega_m( cosmo ), cosmo_omega_l( cosmo ), cosmo_h( cosmo ) );

        /* init cosmo spline rad -> z */
        spl = spline_init( SPLINE_TYPE_CUBIC );
        {
            int i;
            double z, dc, dz;
            double zmin = 1e9, zmax = -1e9;

            /* pad redshift bounds by 10% for spline and z-distortions */
            if( conf->zmin >= 0 ) {
                zmin = 0.9 * conf->zmin;
            } else {
                zmin = 1.1 * conf->zmin;
            }
            zmax = 1.1 * conf->zmax;

            /* these are used for an initial trim! */
            rmin = cosmo_dist_co_los( cosmo, zmin );
            rmax = cosmo_dist_co_los( cosmo, zmax );

            fprintf( stderr, "cosmo> Redshift trimmed to: %g - %g\n", conf->zmin, conf->zmax );
            fprintf( stderr, "cosmo> Creating spline for: %g - %g\n", zmin, zmax );
            fprintf( stderr, "cosmo> Trim comoving distance: %g < r < %g\n", rmin, rmax );

            dz = ( zmax - zmin ) / ( NSPLINE - 1 );
            for( z = zmin, i = 0; i < NSPLINE; i++, z += dz ) {
                dc = cosmo_dist_co_los( cosmo, z );
                spline_data_add( spl, dc, z );  /* it's actually a reverse spline */
            }
            spline_data_finalize( spl );
        }

        /* initialize MANGLE polygon reading */
        if( conf->file_mask != NULL ) {
            fprintf( stderr, "sky> Using mangle polygon file: %s\n", conf->file_mask );
            ply = mply_read_file( conf->file_mask );
            fprintf( stderr, "sky> Downsampling by sky completeness in mask: %s\n",
                     conf->downsample_sky ? "ON" : "OFF" );
            fprintf( stderr, "sky> Minimum weight: %g\n", conf->min_sky_weight );
        } else {
            fprintf( stderr, "sky> No mask, so not by sky completeness in mask.\n" );
            ply = NULL;
        }

        fprintf( stderr, "zsel> Downsampling redshift" );
        if( conf->file_zsel != NULL ) {
            zsel = zsel_read_file( conf->file_zsel );
            fprintf( stderr, " from file: %s\n", conf->file_zsel );
        } else {
            zsel = NULL;
            fprintf( stderr, ": OFF\n" );
        }
    }

    /* Read from ASCII text file, one line at a time */
    sr = sr_init( file_mock );
    fdout = check_fopen( file_out, "w" );
    fprintf( stderr, "\nReading and processing one-by-one ...\n" );
    fprintf( stderr, "  <- %s\n", sr_filename( sr ) );
    fprintf( stderr, "  -> %s\n", file_out );
    while( sr_readline( sr ) ) {
        int i, k, check;
        double x[3], v[3];      /* original coordinates */
        double rx[3], rv[3];    /* remapped coordinates */
        double z, rad, vel;     /* radius and velocity */
        double ran, prob;
        double ra, dec, weight = 1.0;

        if( sr_line_isempty( sr ) )
            continue;

        line = sr_line( sr );
        if( '#' == line[0] )
            continue;

        check = sscanf( line, "%lf %lf %lf %lf %lf %lf", &x[0], &x[1], &x[2], &v[0], &v[1], &v[2] );

        if( check != 6 ) {
            if( check < 3 ) {
                fprintf( stderr,
                         "Error: mock input error on line %d in file: %s\n",
                         sr_linenum( sr ), sr_filename( sr ) );
                exit( EXIT_FAILURE );
            }
            v[0] = v[1] = v[2] = 0.0;
        }

        nread += 1;

        /* Normalize coordinates to [0,1] */
        for( i = 0; i < 3; i++ ) {
            x[i] /= conf->lbox;
            if( x[i] < 0 || x[i] > 1 ) {
                fprintf( stderr,
                         "Error: input out of expected range [0,LBOX = %g] on line %d in file: %s\n",
                         conf->lbox, sr_linenum( sr ), sr_filename( sr ) );
            }
        }

        /* check if we want to flip any axes */
        for( i = 0; i < 3; i++ ) {
            if( conf->pre_rot[i] != 0 ) {
                double ang = conf->pre_rot[i] * ( PI / 2.0 );

                /* center on the origin (0,0,0) */
                for( k = 0; k < 3; k++ )
                    x[k] -= 0.5;
                ct_rotate_about_axis( ang, x, i );
                ct_rotate_about_axis( ang, v, i );
                /* correct recentering to range [0,1] */
                for( k = 0; k < 3; k++ )
                    x[k] += 0.5;
            }
        }

        /* remap point into cuboid, include velocity */
        R.Transform( x[0], x[1], x[2], rx[0], rx[1], rx[2] );
        R.TransformVelocity( v[0], v[1], v[2], rv[0], rv[1], rv[2] );

        /* translate, and recast back to input dimensions */
        for( i = 0; i < 3; i++ ) {
            x[i] = rx[i] * conf->lbox + conf->t[i];
            v[i] = rv[i];
        }

        /* get distance and peculiar velocity of point, these should be rotationally invariant */
        rad = ct_vec_mag( x );
        vel = ct_vec_dot( v, x ) / rad;

        /* let's see if we pass redshift tests */

        /* radial comoving distance -> redshift  */
        if( rad < rmin )
            continue;
        if( rad > rmax )
            continue;

        z = spline_eval( spl, rad );

        /*  add redshift distortion  */
        if( conf->zspace == 1 ) {
            /* need velocity in physical units */
            double dz = vel * ( 1.0 + conf->zin ) / SPEED_OF_LIGHT;
            z += dz;
        }

        /* trim by redshift */
        if( z > conf->zmax )
            continue;
        if( z < conf->zmin )
            continue;

        /* check redshift selection */
        if( zsel != NULL ) {
            prob = zsel_eval( zsel, z );
            if( prob < 1.0 ) {
                ran = rng_uniform( rng );
                if( ran > prob )
                    continue;
            }
        }

        /* so all redshift / radial selection is now done, check sky next */

        /* rotate
         * NOTE: since we added z-distortion already, we don't bother updating
         * velocities by these rotations. DON'T USE VELOCITIES PAST HERE!
         */
        for( i = 0; i < 3; i++ ) {
            if( conf->rot[i] != 0 )
                ct_rotate_about_axis( ct_rad_from_deg( conf->rot[i] ), x, i );
        }

        /* project onto sky: (x,y,z) -> (ra,dec) */
        ct_xyz_to_radec( x, &ra, &dec );

        /* check mangle mask  */
        if( ply != NULL ) {
            MANGLE_INT ipoly;
            ipoly = mply_find_polyindex_radec( ply, ra, dec );

            /* -1 means not in mask.ply */
            if( ipoly < 0 )
                continue;

            /* since we exist, let's look up the associated weight */
            weight = mply_weight_from_index( ply, ipoly );
            if( weight < conf->min_sky_weight )
                continue;

            /* downsample to the weight */
            if( conf->downsample_sky ) {
                ran = rng_uniform( rng );
                if( ran > weight )
                    continue;
            }
        }

        /* output mock point on sky */
        fprintf( fdout, "%10.6f % 10.6f %10.7f %10.7f\n", ra, dec, z, weight );
        fflush( fdout );
        nout += 1;
    }

    fprintf( stderr, "Finished: %zd -> %zd\n", nread, nout );

    /* cleanup / kill */
    fclose( fdout );
    sr_kill( sr );
    rng_kill( rng );
    spline_kill( spl );
    cosmo_kill( cosmo );

    if( ply != NULL ) {
        mply_kill( ply );
    }
    if( zsel != NULL ) {
        zsel_kill( zsel );
    }

    return EXIT_SUCCESS;
}
