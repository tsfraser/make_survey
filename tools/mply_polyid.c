#include <stdlib.h>
#include <stdio.h>

#include <minimal_mangle.c>
#include <simple_reader.c>

int
main( int argc, char **argv )
{
    MANGLE_PLY *ply;
    simple_reader *sr;

    if( argc < 3 ) {
        printf( "Usage: %s  POLYGON  RA_DEC_FILE > OUTPUT \n", argv[0] );
        return EXIT_FAILURE;
    }
    ply = mply_read_file( argv[1] );
    sr = sr_init( argv[2] );

    while( sr_readline( sr ) ) {
        char *line;
        int check;
        MANGLE_INT ipoly;       /* this is an internal index */
        MANGLE_INT polyid;      /* listed POLYID */
        double ra, dec;

        line = sr_line( sr );

        if( sr_line_isempty( sr ) )
            continue;
        if( '#' == line[0] )
            continue;

        check = sscanf( line, "%lf %lf", &ra, &dec );
        if( 2 != check ) {
            fprintf( stderr, "WARNING: skipped line, couldn't read RA/DEC on line %d in file %s\n",
                     sr_linenum( sr ), sr_filename( sr ) );
            continue;
        }

        ipoly = mply_find_polyindex_radec( ply, ra, dec );
        polyid = mply_polyid_from_index( ply, ipoly );

        fprintf( stdout, "%6zd %s\n", ( ssize_t ) polyid, line );
    }

    ply = mply_kill( ply );
    sr = sr_kill( sr );

    return EXIT_SUCCESS;
}
