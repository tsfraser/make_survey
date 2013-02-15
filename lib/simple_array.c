/* simplified array usage, with resizing and sanity checks
 *
 * Cameron McBride
 * cameron.mcbride AT gmail.com
 * January 2013
 * */
#pragma once
#ifndef SIMPLE_ARRAY_DEFINED
#define SIMPLE_ARRAY_DEFINED 1

#include <stdlib.h>
#include <stddef.h>
#include <assert.h>

#include <check_alloc.c>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define SA_GROW_BY_FACTOR 2

typedef struct {                /* member names chosen to match 'man 3 calloc' */
    size_t count;               /* count: number of elements */
    size_t size;                /* size: memory size of each element */
    void *data;                 /* data: pointer to actual storage array */
} simple_array;

static inline int
sa_check_not_null( simple_array const *const sa )
{
    if( NULL == sa ) {
        fprintf( stderr, "Error: simple_array expected to be NOT NULL!\n" );
        assert( NULL == sa );
    }
    return TRUE;
}

void
sa_alloc( simple_array * const sa, const size_t count, const size_t size )
{
    sa_check_not_null( sa );
    sa->data = check_alloc( count, size );
    sa->count = count;
    sa->size = size;
}

void
sa_clean( simple_array * const sa )
{
    if( NULL != sa ) {
        CHECK_FREE( sa->data );
        sa->count = 0;
        sa->size = 0;
    }
}

simple_array *
sa_init( const size_t count, const size_t size )
{
    simple_array *sa;
    sa = ( simple_array * ) check_alloc( 1, sizeof( simple_array ) );
    sa_alloc( sa, count, size );
    return sa;
}

simple_array *
sa_kill( simple_array * sa )
{
    sa_clean( sa );
    CHECK_FREE( sa );
    return NULL;
}

void *
sa_set_length( simple_array * const sa, const size_t new_count )
{
    sa_check_not_null( sa );
    if( new_count != sa->count ) {
        sa->data = check_realloc( sa->data, new_count, sa->size );
        sa->count = new_count;
    }

    return sa->data;
}

void *
sa_ensure_length( simple_array * const sa, const size_t count )
{
    sa_check_not_null( sa );
    size_t new_count;
    /* a no-op if we already think it's big enough */
    if( sa->count >= count )
        return sa->data;

    /* for efficiency, try to realloc to a FACTOR of new elements */
    new_count = sa->count * SA_GROW_BY_FACTOR;

    return sa_set_length( sa, new_count );
}

static inline int
sa_check_length( simple_array const *const sa, const size_t count )
{
    sa_check_not_null( sa );
    if( sa->count >= count )
        return TRUE;
    else
        return FALSE;
}

static inline void *
sa_data( simple_array const *const sa )
{
    sa_check_not_null( sa );
    return sa->data;
}

static inline size_t
sa_len( simple_array const *const sa )
{
    sa_check_not_null( sa );
    return sa->count;
}

#endif
