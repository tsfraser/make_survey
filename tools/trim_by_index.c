#include <stdlib.h>
#include <stdio.h>

#include <simple_reader.c>

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

int
get_next_index( simple_reader * sri, const int last_index )
{
    int index, check;
    char *iline = NULL;

    while( sr_readline( sri ) ) {
        iline = sr_line( sri );
        if( sr_line_isempty( sri ) )
            continue;
        if( '#' == iline[0] )
            continue;
        break;
    }

    if( sr_eof( sri ) || NULL == iline )
        return -1;

    check = sscanf( iline, "%d", &index );
    if( check != 1 ) {
        fprintf( stderr, "ERROR: could not read index on line %d of %s\n",
                 sr_linenum( sri ), sr_filename( sri ) );
        exit( EXIT_FAILURE );
    }

    if( index <= 0 || index <= last_index ) {
        fprintf( stderr, "ERROR: bad index on line %d of %s\n",
                 sr_linenum( sri ), sr_filename( sri ) );
        fprintf( stderr, "  -> INDEX must be positive non-zero and monotonically increasing\n" );
        exit( EXIT_FAILURE );
    }

    return index;
}

int
main( int argc, char **argv )
{
    /* these could also be declared as "void *" */
    simple_reader *srd, *sri;
    int index = -1;

    int nindex = 0, nindex_used = 0;
    int ndata = 0, nkeep = 0;

    if( argc < 3 ) {
        fprintf( stderr, "Usage: %s  FILE  INDEX  >  TRIMMED_FILE\n", argv[0] );
        fprintf( stderr, "  **NOTE: assumes sorted INDEX (monotonically increasing)\n" );
        return EXIT_FAILURE;
    }

    srd = sr_init( argv[1] );   /* input DATA file */
    sri = sr_init( argv[2] );   /* input INDEX file */
    fprintf( stderr, "TRIMMING:   %s\n", sr_filename( srd ) );
    fprintf( stderr, "FROM INDEX: %s\n", sr_filename( sri ) );

    /* loop over input data file */
    while( sr_readline( srd ) ) {
        char *line;

        line = sr_line( srd );

        if( sr_line_isempty( srd ) )
            continue;
        if( '#' == line[0] )
            continue;

        ndata += 1;
        /* ndata is the appropriate 1-based input index, comments and blank lines skipped */

        while( index < ndata ) {
            index = get_next_index( sri, index );
            if( index < 1 )
                break;
            nindex += 1;
        }

        if( ndata == index ) {
            nkeep += 1;
            fprintf( stdout, "%s\n", line );
        }

        if( sr_eof( sri ) )
            break;
    }

    nindex_used = nindex;

    /* check to see how much longer index is */
    while( !sr_eof( sri ) ) {
        index = get_next_index( sri, index );
        if( index < 1 )
            break;
        nindex += 1;
    }
    fprintf( stderr, "INDEX: used %d of %d entries\n", nindex_used, nindex );

    srd = sr_kill( srd );
    sri = sr_kill( sri );

    fprintf( stderr, "INPUT: kept %d of %d read (could be truncated due to short INDEX)\n",
             nkeep, ndata );

    return EXIT_SUCCESS;
}
