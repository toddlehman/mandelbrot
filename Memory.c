/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Memory.h"


//------------------------------------------------------------------------------
public_function
size_t mem_size(void *memory)
{
  if (memory != NULL)
    return malloc_size(memory);
  else
    return 0;
}


//------------------------------------------------------------------------------
public_function
void *mem_alloc(size_t item_count, size_t item_size)
{
  size_t alloc_size = item_count * item_size;
  void *memory = malloc(alloc_size);
  if (!memory)
    error_exit("malloc() failed");
  memset(memory, ~0, alloc_size);
  return memory;
}


//------------------------------------------------------------------------------
public_function
void *mem_alloc_clear(size_t item_count, size_t item_size)
{
  size_t alloc_size = item_count * item_size;
  void *memory = calloc(1, alloc_size);
  if (!memory)
    error_exit("calloc() failed");
  return memory;
}


//------------------------------------------------------------------------------
public_function
void *mem_realloc(void *memory, size_t item_count, size_t item_size)
{
  memory = realloc(memory, item_count * item_size);
  if (!memory)
    error_exit("realloc() failed");
  return memory;
}


//------------------------------------------------------------------------------
public_function
void *mem_alloc_copy(void *memory, size_t item_count, size_t item_size)
{
  assert(memory != NULL);
  void *new_memory = mem_alloc(item_count, item_size);
  if (new_memory != NULL)
    mem_copy(new_memory, memory, item_count, item_size);
  return new_memory;
}


//------------------------------------------------------------------------------
public_function
void *mem_alloc_clone(void *memory)
{
  assert(memory != NULL);
  return mem_alloc_copy(memory, mem_size(memory), 1);
}


//------------------------------------------------------------------------------
public_function
void mem_clear(void *memory, size_t item_count, size_t item_size)
{
  assert(memory != NULL);
  memset(memory, 0, item_count * item_size);
}


//------------------------------------------------------------------------------
public_function
void mem_copy(void *memory_dest, void *memory_src, size_t item_count, size_t item_size)
{
  assert(memory_dest != NULL);
  assert(memory_src != NULL);
  memcpy(memory_dest, memory_src, item_count * item_size);
}


//------------------------------------------------------------------------------
// Note: This routine should not be called directly, but only via the
// mem_dealloc() macro, which is defined in the header file.
public_function
void mem_dealloc_(void **memory)
{
  assert(memory != NULL);
  
  if (*memory)
  {
    free(*memory);
    *memory = NULL;
  }
}


//------------------------------------------------------------------------------
