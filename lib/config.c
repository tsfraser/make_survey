
typedef struct {

    int u[9];                   /* remap matrix */
    double t[3];                /* translation */
    double rot[3];              /* three rotation angles, in degrees */
    double lbox;                /* size of box */
    double omega_m;             /* omega_m */
    double omega_l;             /* omega_l */
    double hubble;              /* dimensionaless hubble: Ho = 100 km/s/Mpc */

    unsigned int seed;          /* used for downsampling */
    char *file_mask;
    char *file_zsel;

    double min_sky_weight;

    int zspace;                 /* output in redshift space? 0 = NO, 1 = YES */
    int downsample_sky;         /* randomly downsample according to sky weight in mask.ply */

    double zin;                 /* mock simulation input redshift, for velocity conversion */
    double zmin, zmax;          /* redshift bounds */

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

    c->lbox = -1.0;
    c->omega_m = 0.0;
    c->omega_l = 0.0;
    c->hubble = 1.0;

    c->file_mask = NULL;
    c->file_zsel = NULL;
    c->seed = 0;
    c->zin = 1e9;               /* BAD FLAG: unphysical */
    c->zmin = 1e20;
    c->zmax = -1e20;

    return c;
}
