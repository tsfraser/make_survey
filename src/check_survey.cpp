/* check_survey
 * ------------
 *
 * GOAL: Perform opposite transformations of make_survey, i.e. RDZ -> XYZ.
 * The primary purpose of this code is to check and verify input to make_survey,
 * so it should read the same configuration input, but reverse operations.
 *
 * The "tricky" parts of the configuration this code should help set are
 * 1. Sky de-projection (using proper cosmological distance)
 * 2. Rotation
 * 3. Check translation
 * 4. Check remapping
 *
 * Cameron McBride
 * cameron.mcbride@gmail.com
 * March 2013
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>

#include <cuboid.h>

#include <simple_reader.c>
#include <check_fopen.c>

#include <config.c>
#include <rng_gsl.c>
#include <spline_gsl.c>
#include <coord_transforms.c>
#include <cosmology.c>

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
    int i;

    char *file_config = NULL;
    char *file_rdz = NULL;

    size_t nread = 0, nout = 0;

    char *line;
    simple_reader *sr;

    CONFIG *conf;
    void *rng;
    COSMO *cosmo;
    SPLINE *spl;
    double bnd_ra[2], bnd_dec[2], bnd_z[2];
    double min[3], max[3];
    double r[3];
    int randz = 0;

    if( argc != 3 ) {
        fprintf( stderr, "USAGE:\t%s  CONFIG_FILE  RDZ_INPUT (Z input optional)\n", argv[0] );
        exit( EXIT_FAILURE );
    }

    file_config = argv[1];
    file_rdz = argv[2];
    conf = conf_readfile( file_config );

    /* initialize cuboid remapping */
    Cuboid R( conf->u );
    {
        r[0] = R.L1 * conf->lbox;
        r[1] = R.L2 * conf->lbox;
        r[2] = R.L3 * conf->lbox;

        fprintf( stdout, "check_survey> ORIGINAL LBOX = %7.2f\n", conf->lbox );
        fprintf( stdout, "check_survey> REMAPPING L1 -> %7.2f ( %7.5f )\n", r[0], R.L1 );
        fprintf( stdout, "check_survey> REMAPPING L2 -> %7.2f ( %7.5f )\n", r[1], R.L2 );
        fprintf( stdout, "check_survey> REMAPPING L3 -> %7.2f ( %7.5f )\n", r[2], R.L3 );
    }

    /* do more initialization */
    {
        /* initialize cosmology */
        cosmo = cosmo_init(  );
        cosmo_set_omega_m( cosmo, conf->omega_m );
        cosmo_set_omega_l( cosmo, conf->omega_l );
        cosmo_set_h( cosmo, 1.0 );

        fprintf( stdout, "check_survey> COSMOLOGY: (Om, Ol, h) = (%g, %g, %.3f)\n",
                 cosmo_omega_m( cosmo ), cosmo_omega_l( cosmo ), cosmo_h( cosmo ) );

        /* init cosmo spline rad -> z */
        spl = spline_init( SPLINE_TYPE_CUBIC );
        {
            double z, dc, dz;
            double zmin = conf->zmin;
            double zmax = conf->zmax;

            fprintf( stdout, "check_survey> COSMOLOGY: creating spline (z -> Dc) for z = %g - %g\n",
                     zmin, zmax );
            dz = ( zmax - zmin ) / ( NSPLINE - 1 );
            for( z = zmin, i = 0; i < NSPLINE; i++, z += dz ) {
                dc = cosmo_dist_co_los( cosmo, z );
                spline_data_add( spl, z, dc );  /* forward spline this time */
            }
            spline_data_finalize( spl );
        }

        if( conf->file_mask ) {
            fprintf( stdout, "check_survey> IGNORING skymask: %s (not trimming input)\n",
                     conf->file_mask );
        }

        for( i = 0; i < 3; i++ ) {
            min[i] = DBL_MAX;
            max[i] = DBL_MAX * -1.0;
        }

        bnd_z[0] = DBL_MAX;
        bnd_ra[0] = DBL_MAX;
        bnd_dec[0] = DBL_MAX;

        bnd_z[1] = DBL_MAX * -1.0;
        bnd_ra[1] = DBL_MAX * -1.0;
        bnd_dec[1] = DBL_MAX * -1.0;

        rng = rng_init( conf->seed );
    }

    fprintf( stdout, "check_survey>   CONFIG: %s\n", file_config );
    fprintf( stdout, "check_survey>   INPUT:  %s\n", file_rdz );

    /* Read from ASCII text file, one line at a time */
    sr = sr_init( file_rdz );
    fprintf( stdout, "check_survey> PROCESSING line-by-line...\n" );
    while( sr_readline( sr ) ) {
        int check;
        double x[3];
        double rad, ra, dec, z;

        if( sr_line_isempty( sr ) )
            continue;

        line = sr_line( sr );
        if( '#' == line[0] )
            continue;

        check = sscanf( line, "%lf %lf %lf", &ra, &dec, &z );
        if( check < 3 ) {
            z = rng_uniform( rng ) * ( conf->zmax - conf->zmin ) + conf->zmin;
            randz = 1;
            if( check < 2 ) {
                fprintf( stderr,
                         "Error: RDZ input error on line %d in file: %s\n",
                         sr_linenum( sr ), sr_filename( sr ) );
                exit( EXIT_FAILURE );
            }
            if( randz < 0 ) {
                fprintf( stderr,
                         "Error: RDZ input error on line %d in file: %s\n",
                         sr_linenum( sr ), sr_filename( sr ) );
                fprintf( stderr, "Error: both 2 and 3 column input (redshift)!\n" );
                exit( EXIT_FAILURE );
            }
        } else {
            if( randz > 0 ) {
                fprintf( stderr,
                         "Error: RDZ input error on line %d in file: %s\n",
                         sr_linenum( sr ), sr_filename( sr ) );
                fprintf( stderr, "Error: both 2 and 3 column input (redshift)!\n" );
                exit( EXIT_FAILURE );
            }
            randz = -1;
        }

        nread += 1;

        /* trim by redshift */
        if( z < conf->zmin || z > conf->zmax )
            continue;

        if( ra < bnd_ra[0] )
            bnd_ra[0] = ra;
        else if( ra > bnd_ra[1] )
            bnd_ra[1] = ra;

        if( dec < bnd_dec[0] )
            bnd_dec[0] = dec;
        else if( dec > bnd_dec[1] )
            bnd_dec[1] = dec;

        if( z < bnd_z[0] )
            bnd_z[0] = z;
        else if( z > bnd_z[1] )
            bnd_z[1] = z;

        /* obtain comoving distance (radius) from redshift */
        rad = spline_eval( spl, z );

        /* project (ra,dec) into cartesian coordinates */
        ct_radec_to_xyz( ra, dec, x );

        /* reverse tranform -- perform in OPPOSITE order and flip sign of angles */
        for( i = 2; i >= 0; i-- )
            if( conf->rot[i] != 0 )
                ct_rotate_about_axis( ct_rad_from_deg( -1.0 * conf->rot[i] ), x, i );

        /* now scale each coordinate by radial distance! */
        for( i = 0; i < 3; i++ ) {
            x[i] *= rad;
            if( x[i] < min[i] )
                min[i] = x[i];
            else if( x[i] > max[i] )
                max[i] = x[i];
        }

        nout += 1;
    }

    fprintf( stdout, "check_survey> DONE: %zd -> %zd\n", nread, nout );
    if( randz > 0 ) {
        fprintf( stdout,
                 "check_survey> REDSHIFT GENERATED: %g < z < %g (uniformly with seed: %u)\n",
                 conf->zmin, conf->zmax, conf->seed );
    } else {
        fprintf( stdout, "check_survey> TRIM REDSHIFT: %g < z < %g\n", conf->zmin, conf->zmax );
    }

    fprintf( stdout, "check_survey> INPUT DATA (TRIMMED): \n" );
    fprintf( stdout, "check_survey>   %-4s MIN,MAX = %10f,%10f\n", "RA:", bnd_ra[0], bnd_ra[1] );
    fprintf( stdout, "check_survey>   %-4s MIN,MAX = %10f,%10f\n", "DEC:", bnd_dec[0], bnd_dec[1] );
    fprintf( stdout, "check_survey>   %-4s MIN,MAX = %10f,%10f\n", "Z:", bnd_z[0], bnd_z[1] );

    fprintf( stdout, "check_survey> USING: rotate  %g,%g,%g (reversed for this check)\n",
             conf->rot[0], conf->rot[1], conf->rot[2] );

    fprintf( stdout, "check_survey> USING: remap  %d,%d,%d  %d,%d,%d  %d,%d,%d\n",
             conf->u[0], conf->u[1], conf->u[2],
             conf->u[3], conf->u[4], conf->u[5], conf->u[6], conf->u[7], conf->u[8] );
    fprintf( stdout, "check_survey> USING: translate  %g,%g,%g\n",
             conf->t[0], conf->t[1], conf->t[2] );
    {
        double c[3];
        int bounded = 1, within_bnds[3];
        double vol = 1, vol_remap = 1;
        char label[][3] = { "X:", "Y:", "Z:" };

        fprintf( stdout, "check_survey> RESULTS:\n" );
        for( i = 0; i < 3; i++ ) {
            double mean = ( max[i] + min[i] ) / 2.0;
            double range = fabs( max[i] - min[i] );
            double t = conf->t[i];
            within_bnds[i] = 0;

            vol *= range;       /* bounding cubical volume */
            vol_remap *= r[i];

            if( min[i] > t )
                if( max[i] < t + r[i] )
                    within_bnds[i] = 1;

            bounded *= within_bnds[i];
            c[i] = round( mean - r[i] / 2.0 );

            fprintf( stdout,
                     "check_survey>  %3s min / remap (pad): %8.2f / %8.2f (%8.2f) %s\n",
                     label[i], min[i], t, min[i] - t, min[i] > t ? "" : "XXX" );
            fprintf( stdout,
                     "check_survey>  %3s max / remap (pad): %8.2f / %8.2f (%8.2f) %s\n",
                     label[i], max[i], r[i] + t, r[i] + t - max[i],
                     max[i] < t + r[i] ? "" : "XXX" );
            fprintf( stdout,
                     "check_survey>  %3s range / remapped:  %8.2f / %8.2f\n",
                     label[i], range, r[i] );
            fprintf( stdout,
                     "check_survey>  %3s mean  / remapped:  %8.2f / %8.2f\n",
                     label[i], mean, ( t + r[i] / 2.0 ) );
            fprintf( stdout, "check_survey>\n" );

        }

        fprintf( stdout, "check_survey> BOUNDARY CHECK" );
        if( bounded ) {
            fprintf( stdout, " OK! (all axes within bounds of remapped volume)\n" );
        } else {
            for( i = 0; i < 3; i++ ) {
                fprintf( stdout, "  %s%s", label[i], within_bnds[i] ? "PASS" : "FAIL" );
            }
            fprintf( stdout, "  (note the 'XXX' above)\n" );
        }

        fprintf( stdout, "check_survey> INPUT BOUNDING VOLUME: %g\n", vol );
        fprintf( stdout, "check_survey> REMAP BOUNDING VOLUME: %g\n", vol_remap );
        fprintf( stdout, "check_survey> SUGGESTED CENTERING: translate  %g,%g,%g\n",
                 c[0], c[1], c[2] );
    }

    sr_kill( sr );
    rng_kill( rng );
    spline_kill( spl );
    cosmo_kill( cosmo );

    return EXIT_SUCCESS;
}
