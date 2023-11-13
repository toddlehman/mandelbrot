/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//------------------------------------------------------------------------------
// ALLOCATION AND DEALLOCATION

extern_public_function
  size_t mem_size(void *memory);

extern_public_function
  void *mem_alloc(size_t item_count, size_t item_size);

extern_public_function
  void *mem_alloc_clear(size_t item_count, size_t item_size);

extern_public_function
  void mem_clear(void *memory, size_t item_count, size_t item_size);

extern_public_function
  void *mem_realloc(void *memory, size_t item_count, size_t item_size);

extern_public_function
  void *mem_alloc_copy(void *memory, size_t item_count, size_t item_size);

extern_public_function
  void *mem_alloc_clone(void *memory);

extern_public_function
  void mem_copy(void *memory_dest, void *memory_src,
                size_t item_count, size_t item_size);

// A note about the memory deallocation routine:  The actual routine, named
// mem_dealloc_(), should not be called directly, as it is wrapped by a macro
// named mem_dealloc(), which allows the caller to write
//
//    mem_dealloc(&pointer)
//
// without getting a stupid compiler warning about incompatible pointer types
// and without having to explicitly write
//
//    mem_dealloc((void **)&pointer)
//
extern_public_function
  void mem_dealloc_(void **memory);

// This works in most cases, but causes problems with ARC in cases of memory
// pointing to Objective-C objects:
//#define  mem_dealloc(x)  mem_dealloc_((void **)(x))

// This works for all cases:
#define  mem_dealloc(x)  mem_dealloc_((void **)(void *)(x))

