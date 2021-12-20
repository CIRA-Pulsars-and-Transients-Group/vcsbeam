/********************************************************
 *                                                      *
 * Licensed under the Academic Free License version 3.0 *
 *                                                      *
 ********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <mwalib.h>

#include "vcsbeam.h"

read_buffer *vmInitReadBuffer( size_t read_size, size_t copy_size )
/*
 *  <----------------------- buffer_size ------------------------->
 *  <-copy_size-> <---------------- read_size -------------------->
 *  <---------------- read_size --------------------> <-copy_size->
 * |-------------|-----------------------------------|-------------|
 * ^             ^                                   ^
 * buffer        read_ptr                            copy_from_ptr
 * copy_to_ptr
 */
{
    // If read size is zero, silently do nothing
    if (read_size == 0)
        return NULL;

    // Allocate memory for the read_buffer struct itself
    read_buffer *rb = (read_buffer *)malloc( sizeof(read_buffer) );

    // Set the size and allocate the host memory
    rb->buffer_size = read_size + copy_size;
    cudaMallocHost( (void **)&(rb->buffer), rb->buffer_size );
    cudaCheckErrors( "vmInitReadBuffer: cudaMallocHost failed" );

    // Set up the "read" and "copy" pointers
    rb->read_size     = read_size;
    rb->copy_size     = copy_size;
    rb->read_ptr      = (void *)((char *)rb->buffer + copy_size);
    rb->copy_from_ptr = (void *)((char *)rb->buffer + read_size);
    rb->copy_to_ptr   = rb->buffer;

    // Initially, set as unlocked
    rb->locked = false;

    // Return pointer to newly allocated read_buffer
    return rb;
}

void vmFreeReadBuffer( read_buffer *rb )
{
    cudaFreeHost( rb->buffer );
    cudaCheckErrors( "vmFreeReadBuffer: cudaFreeHost failed" );

    free( rb );
}

vm_error vmReadBufferCopyMargin( read_buffer *rb )
{
    // Do nothing, in a variety of cases
    if (rb == NULL)
        return VM_READ_BUFFER_NOT_SET;

    if (rb->buffer == NULL)
        return VM_READ_BUFFER_NOT_SET;

    if (rb->copy_from_ptr == NULL ||
            rb->copy_to_ptr == NULL ||
            rb->copy_size == 0)
        return VM_SUCCESS;

    // Execute copy
    memcpy( rb->copy_to_ptr, rb->copy_from_ptr, rb->copy_size );

    return VM_SUCCESS;
}
