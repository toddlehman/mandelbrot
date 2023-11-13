/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//-----------------------------------------------------------------------------
// FIXED-POINT REAL NUMBER STRUCTURE

#pragma pack(push, 4)

typedef struct
{
  bool    s;      // Sign bit (0 = positive, 1 = negative).
  int     k;      // Count of 32-bit values in array.
  uint32  v[0];   // 32-bit values comprising the number.
}
DeepReal;

#pragma pack(pop)


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_constructor
  DeepReal *dr_create(int bits_deep);

extern_public_constructor
  DeepReal *dr_create_int32(int bits_deep, int32 n);

extern_public_constructor
  DeepReal *dr_create_float64(int bits_deep, float64 x);

extern_public_destructor
  void dr_destroy(DeepReal **p_this);

extern_public_method
  void dr_print_debug(DeepReal *this);

extern_public_method
  void dr_print(DeepReal *this);

extern_public_method
  void dr_assign(DeepReal *this, DeepReal *that);

extern_public_method
  void dr_negate(DeepReal *this);

extern_public_method
  float64 dr_get_float64(DeepReal *this);

extern_public_method
  void dr_set_float64(DeepReal *this, float64 x);

extern_public_method
  int32 dr_get_int32(DeepReal *this);

extern_public_method
  void dr_set_int32(DeepReal *this, int32 n);

extern_public_method
  void dr_mul_two(DeepReal *this);

extern_public_method
  void dr_mul_int32(DeepReal *this, int32 n);

extern_public_method
  void dr_div_two(DeepReal *this);

extern_public_method
  void dr_div_int32(DeepReal *this, int32 n);

extern_public_method
  void dr_add(DeepReal *this, DeepReal *that);

extern_public_method
  void dr_sub(DeepReal *this, DeepReal *that);

extern_public_method
  void dr_mul(DeepReal *this, DeepReal *that);

extern_public_method
  void dr_sqr(DeepReal *this);

extern_public_function
  void dr_test(void);
