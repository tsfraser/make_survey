#pragma once
#ifndef SIMPLE_PARAM_INCLUDED
#define SIMPLE_PARAM_INCLUDED

#include <string.h>
#include <ctype.h>

#include <check_alloc.c>
#include <simple_reader.c>

enum simple_param_type {
    SP_TYPE_CHAR = 10,
    SP_TYPE_SINT = 20,
    SP_TYPE_UINT = 30,
    SP_TYPE_FLOAT = 40,
    SP_TYPE_DOUBLE = 50,
    SP_TYPE_STRING = 100,
    SP_TYPE_ARRAY_SINT = 120,
    SP_TYPE_ARRAY_UINT = 130,
    SP_TYPE_ARRAY_FLOAT = 140,
    SP_TYPE_ARRAY_DOUBLE = 150
};

typedef struct {
    int max_tags;
    int ntags;                  /* number of tags */
    int *ptype;                 /* parameter type */
    int *psize;                 /* parameter size IN ELEMENTS (not bytes */
    char **tag;                 /* tag: parameter name */
    char **dval;                /* default value: NULL means required! */
    void **addr;                /* address */
} simple_param;

simple_param *
sp_init( size_t const max_tags )
{
    simple_param *sp;

    sp = ( simple_param * ) check_alloc( 1, sizeof( simple_param ) );

    sp->max_tags = max_tags;
    sp->ntags = 0;
    sp->ptype = ( int * ) check_alloc( max_tags, sizeof( int ) );
    sp->psize = ( int * ) check_alloc( max_tags, sizeof( int ) );
    sp->tag = ( char ** ) check_alloc( max_tags, sizeof( char * ) );
    sp->addr = ( void ** ) check_alloc( max_tags, sizeof( void * ) );
    sp->dval = ( char ** ) check_alloc( max_tags, sizeof( char * ) );

    return sp;
}

simple_param *
sp_kill( simple_param * sp )
{
    int i;

    for( i = 0; i < sp->ntags; i++ ) {
        CHECK_FREE( sp->tag[i] );
        CHECK_FREE( sp->dval[i] );
    }
    CHECK_FREE( sp->ptype );
    CHECK_FREE( sp->psize );
    CHECK_FREE( sp->tag );
    CHECK_FREE( sp->addr );
    CHECK_FREE( sp->dval );

    CHECK_FREE( sp );
    return sp;
}

void
sp_check_not_null( simple_param const *const sp )
{
    if( NULL == sp ) {
        fprintf( stderr, "ERROR: SIMPLE_PARAM used before allocated!\n" );
        abort(  );
    }
}

void
sp_addtag( simple_param * const sp, char const *const tag, void *const addr,
           int const type, int const len, char const *const default_val )
{
    int i, index;
    size_t taglen;
    sp_check_not_null( sp );

    index = sp->ntags;

    if( index >= sp->max_tags ) {
        fprintf( stderr, "ERROR: SIMPLE_PARAM adding tags beyond MAX allowed!\n" );
        abort(  );
    }

    if( NULL == addr ) {
        fprintf( stderr, "ERROR: SIMPLE_PARAM addtag needs non-NULL address for %s\n", tag );
        abort(  );
    }

    taglen = strlen( tag );
    for( i = 0; i < sp->ntags; i++ ) {
        if( strncmp( tag, sp->tag[i], taglen ) == 0 ) {
            fprintf( stderr, "ERROR: SIMPLE_PARAM addtag defined duplicate tag: %s\n", tag );
            abort(  );
        }
    }
    sp->ntags += 1;
    sp->tag[index] = strdup( tag );
    sp->ptype[index] = type;
    sp->psize[index] = len;
    sp->addr[index] = addr;
    if( default_val != NULL )
        sp->dval[index] = strdup( default_val );
}

void
sp_add_char( simple_param * const sp, char const *const tag, void *const addr,
             char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_CHAR, 1, def );
}

void
sp_add_sint( simple_param * const sp, char const *const tag, void *const addr,
             char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_SINT, 1, def );
}

void
sp_add_uint( simple_param * const sp, char const *const tag, void *const addr,
             char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_UINT, 1, def );
}

void
sp_add_float( simple_param * const sp, char const *const tag, void *const addr,
              char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_FLOAT, 1, def );
}

void
sp_add_double( simple_param * const sp, char const *const tag, void *const addr,
               char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_DOUBLE, 1, def );
}

void
sp_add_string( simple_param * const sp, char const *const tag, void *const addr,
               int const len, char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_STRING, len, def );
}

void
sp_add_array_sint( simple_param * const sp, char const *const tag, void *const addr,
                   int const len, char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_ARRAY_SINT, len, def );
}

void
sp_add_array_uint( simple_param * const sp, char const *const tag, void *const addr,
                   int const len, char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_ARRAY_UINT, len, def );
}

void
sp_add_array_float( simple_param * const sp, char const *const tag, void *const addr,
                    int const len, char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_ARRAY_FLOAT, len, def );
}

void
sp_add_array_double( simple_param * const sp, char const *const tag, void *const addr,
                     int const len, char const *const def )
{
    sp_addtag( sp, tag, addr, SP_TYPE_ARRAY_DOUBLE, len, def );
}

int
sp_check_required( simple_param const *const sp )
{
    int i;
    int ntags_set = 0, error = 0;

    sp_check_not_null( sp );

    for( i = 0; i < sp->ntags; i++ ) {
        if( sp->dval[i] != NULL )
            continue;

        if( sp->tag[i][0] != '\0' ) {
            error = 1;
            fprintf( stderr, "ERROR: missing required parameter '%s'\n", sp->tag[i] );
        }

        ntags_set += 1;
    }

    if( error )
        exit( EXIT_FAILURE );

    return ntags_set;
}

void
sp_print_tags( simple_param const *const sp )
{
    int i;

    sp_check_not_null( sp );

    for( i = 0; i < sp->ntags; i++ ) {
        fprintf( stdout, "SIMPLE_PARAM: '%s' of type %d\n", sp->tag[i], sp->ptype[i] );
    }
}

int
sp_setval( simple_param const *const sp, const int index, char *val )
{
    void *ptr;
    int n, res;
    char *str, *tok, *saveptr;

    ptr = sp->addr[index];
    n = sp->psize[index];

    switch ( sp->ptype[index] ) {
        case SP_TYPE_CHAR:
            res = sscanf( val, "%c", ( char * ) ptr );
            break;
        case SP_TYPE_SINT:
            res = sscanf( val, "%d", ( int * ) ptr );
            break;
        case SP_TYPE_UINT:
            res = sscanf( val, "%u", ( unsigned int * ) ptr );
            break;
        case SP_TYPE_FLOAT:
            res = sscanf( val, "%f", ( float * ) ptr );
            break;
        case SP_TYPE_DOUBLE:
            res = sscanf( val, "%lf", ( double * ) ptr );
            break;
        case SP_TYPE_STRING:
            /* the ptr address must be big enough in match input length! */
            strncpy( ( char * ) ptr, val, n - 1 );
            ( ( char * ) ptr )[n - 1] = '\0';
            res = n = 1;        /* strings are different */
            break;
        case SP_TYPE_ARRAY_SINT:
            {
                int *v;
                v = ( int * ) ptr;
                res = 0;
                for( str = val; res < n; str = NULL ) {
                    tok = strtok_r( str, " \t,;", &saveptr );
                    if( tok == NULL )
                        break;
                    res += sscanf( tok, "%d", &v[res] );
                }
            }
            break;
        case SP_TYPE_ARRAY_UINT:
            {
                unsigned int *v;
                v = ( unsigned int * ) ptr;
                res = 0;
                for( str = val; res < n; str = NULL ) {
                    tok = strtok_r( str, " \t,;", &saveptr );
                    if( tok == NULL )
                        break;
                    res += sscanf( tok, "%u", &v[res] );
                }
            }
            break;
        case SP_TYPE_ARRAY_FLOAT:
            {
                float *v;
                v = ( float * ) ptr;
                res = 0;
                for( str = val; res < n; str = NULL ) {
                    tok = strtok_r( str, " \t,;", &saveptr );
                    if( tok == NULL )
                        break;
                    res += sscanf( tok, "%f", &v[res] );
                }
            }
            break;
        case SP_TYPE_ARRAY_DOUBLE:
            {
                double *v;
                v = ( double * ) ptr;
                res = 0;
                for( str = val; res < n; str = NULL ) {
                    tok = strtok_r( str, " \t,;", &saveptr );
                    if( tok == NULL )
                        break;
                    res += sscanf( tok, "%lf", &v[res] );
                }
            }
            break;
        default:
            fprintf( stderr, "ERROR: SIMPLE_PARAM unknown type for %s\n", sp->tag[index] );
            abort(  );
    }

    if( res != n ) {
        fprintf( stderr, "ERROR: parameter '%s': expected %d elements, got only %d\n",
                 sp->tag[index], n, res );
        exit( EXIT_FAILURE );
    }

    return res;
}

void
sp_set_defaults( simple_param const *const sp )
{
    int i;

    sp_check_not_null( sp );

    for( i = 0; i < sp->ntags; i++ ) {

        if( '\0' == sp->tag[i][0] )     /* this is the "processed" flag */
            continue;
        if( NULL == sp->dval[i] )       /* this is the required flag */
            continue;

        sp_setval( sp, i, sp->dval[i] );
    }
}

int
sp_parse_line( simple_param const *const sp, char *line )
{
    int i, res = 0;
    char *tag;
    size_t ltag, lline;

    sp_check_not_null( sp );

    /* this does not have to be efficient! */
    lline = strlen( line );
    for( i = 0; i < sp->ntags; i++ ) {
        tag = sp->tag[i];
        ltag = strlen( tag );

        if( '\0' == tag[0] )
            continue;

        if( strncmp( line, tag, ltag ) == 0 ) {
            size_t k;
            char *val;

            val = &line[ltag];
            /* skip over initial whitespace in value */
            for( k = ltag; k < lline; k++ ) {
                val = &line[k];
                if( !isblank( ( int ) val[0] ) )
                    break;
            }
            res = sp_setval( sp, i, val );
            tag[0] = '\0';      /* destructive: flag as processed! */
            break;
        }
    }

    return res;
}

int
sp_parse_file( simple_param const *const sp, char *filename )
{
    simple_reader *sr;
    int ntags_read = 0;

    sr = sr_init( filename );
    while( sr_readline( sr ) ) {
        char *line;
        int check;

        line = sr_line( sr );

        if( sr_line_isempty( sr ) )
            continue;
        if( line[0] == '#' )
            continue;

        check = sp_parse_line( sp, line );
        ntags_read += 1;

        if( check < 1 ) {
            fprintf( stderr, "ERROR: unknown or duplicate parameter on line %d in file: %s\n",
                     sr_linenum( sr ), sr_filename( sr ) );
            exit( EXIT_FAILURE );
        }
    }
    sr_kill( sr );

    return ntags_read;
}

#endif
