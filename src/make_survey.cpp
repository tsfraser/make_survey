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
  int ngaltot;
  float nbartot, wfkp, nbar_current;

    char *file_config = NULL;
    char *file_mock = NULL;
    char *file_out = NULL;
    FILE *fdout, *fdi;

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

    fdi = NULL;
    ply = NULL;
    rng = NULL;
    zsel = NULL;

    if( argc != 4 ) {
        fprintf( stderr, "USAGE:\t%s  CONFIG_FILE  MOCK_IN  RDZ_OUT\n", argv[0] );
        exit( EXIT_FAILURE );
    }

    file_config = argv[1];
    file_mock = argv[2];
    file_out = argv[3];
    if( strcmp( file_config, file_out ) == 0 ) {
        fprintf( stderr, "Error: OUTPUT file must be different than INPUT!\n" );
        exit( EXIT_FAILURE );
    }
    if( strcmp( file_mock, file_out ) == 0 ) {
        fprintf( stderr, "Error: OUTPUT file must be different than INPUT!\n" );
        exit( EXIT_FAILURE );
    }

    conf = conf_readfile( file_config );

    fprintf( stderr, "make_survey> PRE-ROTATION about each axis in 90 deg units [%d,%d,%d]\n",
             conf->pre_rot[0], conf->pre_rot[1], conf->pre_rot[2] );

    /* initialize cuboid remapping */
    Cuboid R( conf->u );
    {
        double r[3];
        // print out lattice vectors
        fprintf( stderr, "make_survey> REMAP LATTICE VECTORS:\n" );
        for( int i = 0; i < 3; i++ ) {
            fprintf( stderr, "make_survey>   u%d = (", i + 1 );
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

        fprintf( stdout, "make_survey> ORIGINAL LBOX = %7.2f\n", conf->lbox );
        fprintf( stdout, "make_survey> REMAPPING L1 -> %7.2f ( %7.5f )\n", r[0], R.L1 );
        fprintf( stdout, "make_survey> REMAPPING L2 -> %7.2f ( %7.5f )\n", r[1], R.L2 );
        fprintf( stdout, "make_survey> REMAPPING L3 -> %7.2f ( %7.5f )\n", r[2], R.L3 );

        fprintf( stderr, "make_survey> TRANSLATE X by %8.2f\n", conf->t[0] );
        fprintf( stderr, "make_survey> TRANSLATE Y by %8.2f\n", conf->t[1] );
        fprintf( stderr, "make_survey> TRANSLATE Z by %8.2f\n", conf->t[2] );

        fprintf( stderr, "make_survey> ROTATE about X,Y,Z in degrees: %g,%g,%g\n",
                 conf->rot[0], conf->rot[1], conf->rot[2] );
    }

    /* do more initialization */
    {
        /* initialize cosmology */
        cosmo = cosmo_init(  );
        cosmo_set_omega_m( cosmo, conf->omega_m );
        cosmo_set_omega_l( cosmo, conf->omega_l );
        cosmo_set_h( cosmo, conf->hubble );
        cosmo_set_w0( cosmo, conf->w0);
        cosmo_set_wa( cosmo, conf->wa);

        fprintf( stderr, "make_survey> COSMOLOGY: (Om, Ol, h, w0, wa) = (%g, %g, %.3f,%.3f,%.3f)\n",
                 cosmo_omega_m( cosmo ), cosmo_omega_l( cosmo ), cosmo_h( cosmo ),cosmo_w0( cosmo ), cosmo_wa( cosmo) );

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
            /*printf("rmin,rmax %4.4f %4.4f zmin,zmax %4.4f %4.4f\n",rmin,rmax,zmin,zmax);*/
            fprintf( stderr, "make_survey> COSMOLOGY: creating spline (Dc -> z) for z = %g - %g\n",
                     zmin, zmax );
            fprintf( stderr, "make_survey> REDSHIFT SPACE OUTPUT: %s\n",
                     conf->zspace ? "ON" : "OFF" );
            fprintf( stderr, "make_survey> TRIM REDSHIFT: %g < z < %g\n", conf->zmin, conf->zmax );
            fprintf( stderr, "make_survey> TRIM DISTANCE [trim redshift with a 10 percent buffer!]: %g < r < %g\n", rmin, rmax );

            dz = ( zmax - zmin ) / ( NSPLINE - 1 );
            for( z = zmin, i = 0; i < NSPLINE; i++, z += dz ) {
                dc = cosmo_dist_co_los( cosmo, z );
                /*printf("Spline: dc and z %4.4lf %4.4f\n",dc,z);*/
                spline_data_add( spl, dc, z );  /* it's actually a reverse spline */
            }
            spline_data_finalize( spl );
        }

        /* initialize MANGLE polygon reading */
        if( NULL == conf->file_mask || '\0' == conf->file_mask[0] ) {
            fprintf( stderr, "make_survey> NO SKY MASK: ignoring minimum completeness\n" );
            conf->downsample_sky = OFF;
            ply = NULL;
        } else {
            fprintf( stderr, "make_survey> TRIM SKY: remove regions with completeness < %g\n",
                     conf->min_sky_weight );
            fprintf( stderr, "make_survey> MANGLE POLYGON MASK: %s\n", conf->file_mask );
            ply = mply_read_file( conf->file_mask );
            fprintf( stderr, "make_survey> DOWNSAMPLING by sky completeness in mask: %s\n",
                     conf->downsample_sky ? "ON" : "OFF" );
        }

        fprintf( stderr, "make_survey> DOWNSAMPLING redshift" );
        if( NULL == conf->file_zsel || '\0' == conf->file_zsel[0] ) {
            zsel = NULL;
            fprintf( stderr, ": OFF\n" );
        } else {
            zsel = zsel_read_file( conf->file_zsel );
            fprintf( stderr, " from file: %s\n", conf->file_zsel );
        }

        if( conf->downsample_sky || zsel != NULL ) {
            fprintf( stderr, "make_survey> RANDOM NUMBER GENERATOR initialized with seed = %u\n",
                     conf->seed );
            rng = rng_init( conf->seed );
        }

    }

    

    sr = sr_init( file_mock );
    fdout = check_fopen( file_out, "w" );
    fprintf( stderr, "make_survey>   CONFIG: %s\n", file_config );
    fprintf( stderr, "make_survey>   INPUT:  %s\n", sr_filename( sr ) );
    fprintf( stderr, "make_survey>   OUTPUT: %s\n", file_out );

    // JLT - how many lines in the box?
    ngaltot = 0;
    while(sr_readline(sr))ngaltot++;
    nbartot = ngaltot/pow(conf->lbox,3.0);
    fprintf(stderr,"make_survey> total mock gals: %d nbar=%e\n",ngaltot,nbartot);
    sr = sr_init( file_mock );

    if( conf->make_info ) {
        char *file_info;
        size_t len;

        len = strlen( file_out );
        file_info = ( char * ) check_alloc( len + 10, sizeof( char ) );
        sprintf( file_info, "%s.info", file_out );
        fdi = check_fopen( file_info, "w" );
        fprintf( fdi, "# INDEX FROM: %s (1-based)\n", file_mock );
        fprintf( fdi, "# FOR OUTPUT: %s\n", file_out );
        fprintf( fdi, "# USING CONF: %s\n", file_config );
        fprintf( stderr, "make_survey>   INFO:   %s\n", file_info );
        CHECK_FREE( file_info );
    }

    /* Read from ASCII text file, one line at a time */
    fprintf( stderr, "make_survey> PROCESSING line-by-line...\n" );
    while( sr_readline( sr ) ) {
        int i, k, ncols_in;
        double x[3], v[3];      /* original coordinates */
        double rx[3], rv[3];    /* remapped coordinates */
        double z, rad, vel;     /* radius and velocity */
        double delta_z, z_real, z_red;
        double ra, dec, weight = 1.0;
        double ran1 = 0.0, ran2 = 0.0, prob = 1.0;
        int flag_sat, id_halo;
        double mass_halo;

        if( sr_line_isempty( sr ) )
            continue;

        line = sr_line( sr );
        if( '#' == line[0] )
            continue;

        ncols_in = sscanf( line, "%lf %lf %lf %lf %lf %lf %lf %i %i",
                           &x[0], &x[1], &x[2], &v[0], &v[1], &v[2],
                           &mass_halo, &flag_sat, &id_halo );

        if( ncols_in < 6 ) {
            if( ncols_in < 3 ) {
                fprintf( stderr,
                         "Error: mock input error on line %d in file: %s\n",
                         sr_linenum( sr ), sr_filename( sr ) );
                exit( EXIT_FAILURE );
            }
            v[0] = v[1] = v[2] = 0.0;
        }

        nread += 1;

        /* make RNGs stable: same input catalog and RNG SEED will 
         * now produce same RNG for each object */
        if( zsel != NULL )
            ran1 = rng_uniform( rng );

        if( conf->downsample_sky )
            ran2 = rng_uniform( rng );

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
        /*printf("N_Read: %4.2lu \n",nread);
        printf("XYZ %4.4f %4.4f %4.4f \n",x[0],x[1],x[2]);
        printf("get dist and pecvel floats: %4.2f %4.2f \n", rad,vel);*/
        ct_xyz_to_radec(x,&ra,&dec);
        
        /* let's see if we pass redshift tests */

        /* radial comoving distance -> redshift  */
        if( rad < rmin )
            continue;
        if( rad > rmax )
            continue;

        z_real = spline_eval( spl, rad );
        /*printf ("get z_real float: %4.4f \n",z_real);*/
        /* redshift distortion: need velocity in physical units */
        delta_z = vel * ( 1.0 + conf->zin ) / SPEED_OF_LIGHT;
        z_red = z_real + delta_z;
        /*printf("rsd redshift:%4.4f \n",z_red);*/
        if( conf->zspace )
            z = z_red;
        else
            z = z_real;
        /*printf("just trying radecz from coords... %4.4f %4.4f %4.4f \n",ra,dec,z);*/
        /* trim by redshift */
        if( z > conf->zmax )
            continue;
        if( z < conf->zmin )
            continue;

        /* check redshift selection */
        if( zsel != NULL ) {
            prob = zsel_eval( zsel, z );
            if( prob < 1.0 && ran1 > prob ) {
                continue;
            }
        }
	nbar_current = nbartot*prob;

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
        /*printf("RA,DEC,Z: %4.2f %4.2f %4.4f \n",ra, dec, z);*/
        /* check mangle mask  */
        if( ply != NULL ) {
            MANGLE_INT ipoly;
            ipoly = mply_find_polyindex_radec( ply, ra, dec );

            /* -1 means not in mask.ply */
            if( ipoly < 0 )
                continue;

            /* since we exist, let's look up the associated completeness weight */
            weight = mply_weight_from_index( ply, ipoly );
            if( weight < conf->min_sky_weight )
                continue;

            /* downsample to the weight */
            if( conf->downsample_sky ) {
                if( ran2 > weight )
                    continue;
            }
        }

        /* output mock point on sky */
        //fprintf( fdout, "%10.6f % 10.6f %10.7f\n", ra, dec, z );
	// JLT adds stuff
	wfkp = 1/(1+nbar_current*conf->powspec);


        fprintf( fdout, "%10.6f %10.6f %10.7f %10.7f  %.5e  %.5e\n", ra, dec, z ,weight, nbar_current, wfkp);
        if( conf->make_info > 0 ) {
            if( conf->make_info >= 1 ) {
                fprintf( fdi, "%10zu", nread ); /* this doubles as index */
            }
            if( conf->make_info >= 2 ) {
                fprintf( fdi, " %10.7f", weight );
            }

            /* the additional columns only make sense if they exist in input! */
            if( conf->make_info >= 4 && ncols_in > 3 ) {
                /* calculations require input velocity */
                fprintf( fdi, " %10.7f", z_real );
                fprintf( fdi, " %10.7f", z_red );
            }
            if( conf->make_info >= 6 && ncols_in > 7 ) {
                /* info directly from input file */
                fprintf( fdi, " %12.6e", mass_halo );
                fprintf( fdi, " %1i", flag_sat );
            }
            if( conf->make_info >= 7 && ncols_in > 8 ) {
                /* info directly from input file */
                fprintf( fdi, " %10i", id_halo );
            }
            fprintf( fdi, "\n" );
        }

        nout += 1;
    }

    fprintf( stderr, "make_survey> DONE: %zd -> %zd\n", nread, nout );

    /* cleanup / kill */
    fclose( fdout );
    if( conf->make_info )
        fclose( fdi );
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
