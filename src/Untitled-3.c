
#include <check_alloc.c>
#include <simple_param.c>

typedef struct {

    int u[9];                   /* remap matrix */
    double t[3];                /* translation */
    double rot[3];              /* three rotation angles, in degrees */
    int pre_rot[3];             /* pre-remap axis flips: integral number of 90 deg rotations
                                 * to stare down different faces */

    unsigned int seed;          /* used for downsampling */

    double lbox;                /* size of box */
    double omega_m;             /* omega_m */
    double omega_l;             /* omega_l */
    double w0;                  /* w0 from LCDM */
    double wa;                  /* wa from LCDM */
    double hubble;              /* dimensionaless hubble: Ho = 100 km/s/Mpc */
    double zin;                 /* mock simulation input redshift, for velocity conversion */

    double zmin, zmax;          /* redshift bounds */
    int zspace;                 /* output in redshift space? 0 = NO, 1 = YES */

    char *file_mask;
    int downsample_sky;         /* randomly downsample according to sky weight in mask.ply */
    double min_sky_weight;      /* minimum sky weight to keep */

    char *file_zsel;            /* redshift selection file */

  double powspec;  /* power spectrum value for the FKP weights */

    int make_info;              /* Include additional "Value-Added Columns" in .info file
                                   index:      line number of input mock (1-based)
                                   comp:       completeness (weight in MANGLE file), range is 0-1 typically
                                   z_real:     real redshift (no distortion)
                                   z_redshift: redshift with z-distortion
                                   mass_halo:  mass of parent central halo in Msun units
                                   flag_sat:   is a satellite galaxy?  (0 = no, 1 = yes)
                                   id_halo:    ID of parent halo from input halo catalog
                                   make_info = 0 : DO NOT CREATE
                                   make_info = 1 : index only
                                   make_info = 2 : index, comp
                                   make_info = 4 : index, comp, z_real, z_redshift
                                   make_info = 6 : index, comp, z_real, z_redshift, mass_halo, flag_sat
                                   make_info = 7 : index, comp, z_real, z_redshift, mass_halo, flag_sat, id_halo
                                 */

} CONFIG;

CONFIG *
conf_init( void )
{
    int i, k;
    CONFIG *c;
    c = ( CONFIG * ) check_alloc( 1, sizeof( CONFIG ) );

    for( i = 0; i < 3; i++ ) {
        c->t[i] = 0.0;
        c->rot[i] = 0.0;
        for( k = 0; k < 3; k++ ) {
            if( i == k )
                c->u[i * 3 + k] = 1;
            else
                c->u[i * 3 + k] = 0;
        }
    }

    c->rot[0] = 0.0;
    c->rot[1] = 0.0;
    c->rot[2] = 0.0;

    c->pre_rot[0] = 0;
    c->pre_rot[1] = 0;
    c->pre_rot[2] = 0;

    c->lbox = -1.0;
    c->omega_m = 0.0;
    c->omega_l = 0.0;
    c->w0 = 0.0;
    c->wa = 0.0;
    c->hubble = 1.0;
    c->powspec = 0.0;

    c->file_mask = NULL;
    c->file_zsel = NULL;
    c->seed = 0;
    c->zin = 1e9;               /* BAD FLAG: unphysical */
    c->zmin = 1e20;
    c->zmax = -1e20;

    return c;
}

CONFIG *
conf_kill( CONFIG * c )
{
    CHECK_FREE( c->file_mask );
    CHECK_FREE( c->file_zsel );
    CHECK_FREE( c );
    return c;
}

CONFIG *
conf_readfile( char *filename )
{
    CONFIG *c;

    int n;
    simple_param *sp;
    const size_t filename_length = 1024;

    /* using this string parsing is a bad hack */
    c = conf_init(  );
    c->file_mask = ( char * ) check_alloc( filename_length, sizeof( char ) );
    c->file_zsel = ( char * ) check_alloc( filename_length, sizeof( char ) );

    sp = sp_init( 19 );

    sp_add_array_sint( sp, "remap", c->u, 9, "1 0 0  0 1 0  0 0 1" );
    sp_add_array_sint( sp, "pre_rotate", c->pre_rot, 3, "0 0 0" );

    sp_add_array_double( sp, "translate", c->t, 3, "0 0 0" );
    sp_add_array_double( sp, "rotate", c->rot, 3, "0 0 0" );

    sp_add_string( sp, "file_skymask", c->file_mask, filename_length, "\0" );
    sp_add_string( sp, "file_redshift_selection", c->file_zsel, filename_length, "\0" );

    sp_add_double( sp, "lbox", &c->lbox, NULL );
    sp_add_double( sp, "omega_m", &c->omega_m, NULL );
    sp_add_double( sp, "omega_l", &c->omega_l, NULL );
    sp_add_double( sp, "w0", &c->w0, NULL);
    sp_add_double( sp, "wa", &c->wa, NULL);
    sp_add_double( sp, "hubble", &c->hubble, NULL );
    sp_add_double( sp, "redshift_input", &c->zin, NULL );
    sp_add_double( sp, "max_redshift", &c->zmax, NULL );
    sp_add_double( sp, "min_redshift", &c->zmin, "0.0" );
    sp_add_double( sp, "min_sky_weight", &c->min_sky_weight, "0.0" );

    //JLT
    sp_add_double( sp, "powspec", &c->powspec, NULL );

    sp_add_sint( sp, "downsample_sky", &c->downsample_sky, "0" );
    sp_add_uint( sp, "random_seed", &c->seed, "-12345" );
    sp_add_sint( sp, "redshift_space", &c->zspace, NULL );
    sp_add_sint( sp, "make_info", &c->make_info, "0" );

    sp_set_defaults( sp );
    n = sp_parse_file( sp, filename );
    sp_check_required( sp );

    sp_kill( sp );

    return c;
}
