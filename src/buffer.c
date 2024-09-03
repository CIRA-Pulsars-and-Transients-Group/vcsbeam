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
#include "gpu_macros.h"

/**
 * Initialise a buffer for reading in voltage data.
 *
 * @param read_size The size of the "read" portion of the buffer into which
 *                  the data are read
 * @param copy_size The size of the "copy" portion of the buffer
 *
 * @return A pointer to a newly allocated read buffer
 *
 * The read buffer is designed for situations where the data must be read in
 * a section at a time, but where operations on the data require it to be
 * treated list a continuous stream. In this regard, it is similar to a ring
 * buffer. Here, however, the problem is solved by copying a portion of data
 * at the end of the buffer to the beginning of the buffer at each successive
 * read. For this to be efficient, the size of the needed overlap between
 * consecutive sections must be small enough that copying the data around in
 * this way does not carry a significant amount of overhead.
 *
 * The relationship between the total buffer, the read buffer, and the copy
 * buffer are illustrated in the following diagram.
 *
 * ```
 *  <----------------------- buffer_size ------------------------->
 *                <---------------- read_size -------------------->
 *  <-copy_size->                                     <-copy_size->
 *
 * |-------------|-----------------------------------|-------------|
 * ^             ^                                   ^
 * buffer        read_ptr                            copy_from_ptr
 * copy_to_ptr
 * ```
 *
 * \see vmFreeReadBuffer()
 * \see vmReadBufferCopyMargin()
 * \see vmReadNextSecond()
 */
host_buffer *vmInitReadBuffer( size_t read_size, size_t copy_size )
{
    // If read size is zero, silently do nothing
    if (read_size == 0)
        return NULL;

    // Allocate memory for the host_buffer struct itself
    host_buffer *rb = (host_buffer *)malloc( sizeof(host_buffer) );

    // Set the size and allocate the host memory
    rb->buffer_size = read_size + copy_size;
    gpuMallocHost( (void **)&(rb->buffer), rb->buffer_size );

    // Set up the "read" and "copy" pointers
    rb->read_size     = read_size;
    rb->copy_size     = copy_size;
    rb->read_ptr      = (void *)((char *)rb->buffer + copy_size);
    rb->copy_from_ptr = (void *)((char *)rb->buffer + read_size);
    rb->copy_to_ptr   = rb->buffer;

    // Initially, set as unlocked
    rb->locked = false;

    // Return pointer to newly allocated host_buffer
    return rb;
}

/**
 * Frees the memory allocated in vmInitReadBuffer()
 *
 * @param rb The read buffer to be freed
 *
 * \see vmInitReadBuffer()
 */
void vmFreeReadBuffer( host_buffer *rb )
{
    gpuHostFree( rb->buffer );

    free( rb );
}

/**
 * Copy the data from the end of the read buffer to the beginning.
 *
 * @param rb The read buffer to be operated on
 * @return VM_READ_BUFFER_NOT_SET if either `rb` or `rb&rarr;buffer` is NULL,
 *         or VM_SUCCESS if the memory is copied successfully.
 *
 * This copies `copy_size` bytes from `copy_from_ptr` to `copy_to_ptr`.
 * See vmInitReadBuffer() for a relevant diagram.
 *
 * \see vmInitReadBuffer()
 * \see vmReadNextSecond()
 */
vm_error vmReadBufferCopyMargin( host_buffer *rb )
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
