/* a simple ascii / text file reader interface
 *
 * Cameron McBride
 * cameron.mcbride AT gmail.com
 * January 2013
 * */
#pragma once
#ifndef SIMPLE_READER_DEFINED
#define SIMPLE_READER_DEFINED 1

#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#include <check_fopen.c>
#include <check_alloc.c>

typedef struct {
    FILE *fp;
    char *filename;
    char *buf;
    size_t buf_size;
    size_t buf_len;
    size_t line_num;
} simple_reader;

static inline int
sr_check_not_null( simple_reader const *const sr )
{
    if( NULL == sr ) {
        fprintf( stderr, "Error: simple_reader expected to be NOT NULL!\n" );
        assert( NULL == sr );
    }
    return TRUE;
}

void
sr_buf_alloc( simple_reader * sr, const size_t maxlen )
{
    sr_check_not_null( sr );
    sr->buf_size = maxlen + 2;  /* one for \0 and one to replace possible "\n" */
    sr->buf = ( char * ) check_realloc( sr->buf, sr->buf_size, sizeof( char ) );
}

void
sr_buf_clean( simple_reader * sr )
{
    if( NULL != sr ) {
        CHECK_FREE( sr->buf );
        sr->buf_size = 0;
    }
}

static inline void
sr_set_line_maxlen( simple_reader * sr, const size_t maxlen )
{
    /* just a more sane external name for same functionality */
    sr_buf_alloc( sr, maxlen );
}

size_t
sr_buf_maxlen( simple_reader const *const sr )
{
    sr_check_not_null( sr );
    if( sr->buf_size >= 2 )
        return sr->buf_size - 2;
    return 0;
}

void
sr_buf_wipe( simple_reader * sr )
{
    size_t i;
    for( i = 0; i < sr->buf_size; i++ ) {
        sr->buf[i] = '\0';
    }
}

void
sr_setfile( simple_reader * const sr, char const *const filename )
{
    size_t len;

    if( sr->fp != NULL )
        fclose( sr->fp );
    sr->fp = check_fopen( filename, "r" );

    /* if sr->filename is NULL then check_realloc will calloc it */
    len = strlen( filename ) + 1;
    sr->filename = ( char * ) check_realloc( sr->filename, len, sizeof( char ) );
    strncpy( sr->filename, filename, len );

    sr->line_num = 0;
    sr_buf_wipe( sr );
}

void
sr_alloc( simple_reader * const sr, char const *const filename )
{
    sr_check_not_null( sr );

    if( filename != NULL ) {
        sr_setfile( sr, filename );
    }

    /* choose a reasonable default line buffer size if not already set */
    if( NULL == sr->buf || sr->buf_size < 1 )
        sr_buf_alloc( sr, 1024 );
}

void
sr_clean( simple_reader * const sr )
{
    if( NULL != sr ) {
        fclose( sr->fp );
        CHECK_FREE( sr->filename );
        sr_buf_clean( sr );
        sr->line_num = 0;
    }
}

simple_reader *
sr_init( char const *const filename )
{
    simple_reader *sr;
    sr = ( simple_reader * ) check_alloc( 1, sizeof( simple_reader ) );
    sr_alloc( sr, filename );
    return sr;
}

simple_reader *
sr_kill( simple_reader * sr )
{
    sr_clean( sr );
    CHECK_FREE( sr );
    return NULL;
}

static inline char *
sr_line( simple_reader const *const sr )
{
    sr_check_not_null( sr );
    return sr->buf;
}

static inline int
sr_line_isempty( simple_reader const *const sr )
{
    sr_check_not_null( sr );

    if( NULL == sr->buf )
        return TRUE;
    if( '\0' == sr->buf[0] )
        return TRUE;

    assert( sr->buf_len > 0 );
    return FALSE;
}

static inline int
sr_linenum( simple_reader const *const sr )
{
    sr_check_not_null( sr );
    return sr->line_num;
}

static inline size_t
sr_linelen( simple_reader const *const sr )
{
    sr_check_not_null( sr );
    return sr->buf_len;
}

static inline char *
sr_filename( simple_reader const *const sr )
{
    sr_check_not_null( sr );
    return sr->filename;
}

int
sr_eof( simple_reader const *const sr )
{
    sr_check_not_null( sr );
    return feof( sr->fp );
}

char *
sr_readline( simple_reader * const sr )
{
    char *res;
    size_t l = 0;

    sr_check_not_null( sr );

    if( sr->buf_size < 2 ) {
        fprintf( stderr, "Error: buffer too small! (max_length = %zd)\n", sr->buf_size - 2 );
        exit( EXIT_FAILURE );
    }

    res = fgets( sr->buf, sr->buf_size, sr->fp );
    if( NULL == res ) {         /* check for error or EOF */
        if( feof( sr->fp ) )
            return NULL;

        fprintf( stderr, "Error: Cannot read beyond line %zd of file: %s\n",
                 sr->line_num, sr->filename );

        if( ferror( sr->fp ) )
            perror( "Error:" );

        exit( EXIT_FAILURE );
    }

    sr->line_num += 1;

    l = strcspn( sr->buf, "\n" );
    if( sr->buf_size <= l + 1 ) {
        fprintf( stderr, "Error: exceeded maximum length (%zd) for line %zd in file: %s\n",
                 sr_buf_maxlen( sr ), sr->line_num, sr->filename );
        exit( EXIT_FAILURE );
    }
    sr->buf[l] = '\0';
    sr->buf_len = l;

    return sr->buf;
}

#endif
