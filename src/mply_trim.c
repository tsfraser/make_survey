#include <stdlib.h>
#include <stdio.h>

#include <minimal_mangle.c>
#include <simple_reader.c>

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

int
main( int argc, char **argv )
{
    /* these could also be declared as "void *" */
    MANGLE_PLY *ply;
    simple_reader *sr;

    int reverse_trim = FALSE;
    double min_weight = 0.0;
    size_t nread = 0, nkeep = 0;

    if( argc < 3 ) {
        printf( "Usage: %s  RA_DEC_FILE POLYGON  [MIN_WEIGHT]  [REVERSE_TRIM]  >  OUTPUT\n",
                argv[0] );
        return EXIT_FAILURE;
    }

    fprintf( stderr, "READING polygon file: %s\n", argv[2] );
    ply = mply_read_file( argv[2] );

    if( argc > 3 ) {
        min_weight = strtod( argv[3], NULL );
    }
    if( argc > 4 ) {
        char c, *s;
        s = argv[4];
        c = s[0];

        /* start with it's numeric value */
        sscanf( s, "%d", &reverse_trim );

        if( c == 'y' || c == 'Y' || c == 't' || c == 'T' )
            reverse_trim = TRUE;
        else if( strncasecmp( "on", s, 2 ) == 0 )
            reverse_trim = TRUE;
    }

    if( reverse_trim )
        fprintf( stderr, "FILTERING: vetoing weight >= %g (REVERSED!)\n", min_weight );
    else
        fprintf( stderr, "FILTERING: keeping weight >= %g\n", min_weight );

    sr = sr_init( argv[1] );    /* simple line-by-line reader */
    fprintf( stderr, "PROCESSING: ra dec from %s\n", sr_filename( sr ) );
    while( sr_readline( sr ) ) {
        char *line;
        int check;
        MANGLE_INT index;
        double ra, dec, weight;

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

        nread += 1;

        index = mply_find_polyindex_radec( ply, ra, dec );

        if( index < 0 )
            weight = 0.0;
        else
            weight = mply_weight_from_index( ply, index );

        {
            int skip = FALSE;

            if( weight < min_weight )
                skip = TRUE;

            if( reverse_trim )
                skip = skip ? FALSE : TRUE;

            if( skip )
                continue;
        }

        nkeep += 1;
        fprintf( stdout, "%s\n", line );
    }

    ply = mply_kill( ply );
    sr = sr_kill( sr );

    fprintf( stderr, "DONE: %zu -> %zu\n", nread, nkeep );

    return EXIT_SUCCESS;
}
