/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

// Standard library
#import <stdbool.h>
#import <stdint.h>
#import <stdlib.h>
#import <stdio.h>
#import <inttypes.h>
#import <unistd.h>
#import <time.h>
#import <memory.h>
#import <malloc/malloc.h>
#import <string.h>
#import <ctype.h>
#import <math.h>
#import <complex.h>
#import <limits.h>
#import <assert.h>


//------------------------------------------------------------------------------
// The overloading of the word "static" in C is annoying.  Sometimes it actually
// means static, but many times it means private, and sometimes it is even
// needed for public functions (namely, when declaring inline functions in a
// header file).  For easier reading, the following are used in this program
// instead.

#define  private_static             static
#define  public_static

#define  private_const              static const
#define  public_const               static const

#define  private_function           static
#define  public_function
#define  extern_public_function     extern

#define  private_constructor        static
#define  public_constructor
#define  extern_public_constructor  extern

#define  private_destructor         static
#define  public_destructor
#define  extern_public_destructor   extern

#define  private_method             static
#define  public_method  
#define  extern_public_method       extern

#define  private_inline_function    static inline
#define  public_inline_function     static inline

#define  private_inline_method      static inline
#define  public_inline_method       static inline


//-----------------------------------------------------------------------------
// USEFUL MATH CONSTANTS

#define  PI   M_PI
#define  TAU  (M_PI * 2)


//-----------------------------------------------------------------------------
// USEFUL MACROS

#define  unless(x)  if (!(x))

#define  until(x)   while (!(x))

#define  ELEMENT_COUNT(x)   (sizeof(x) / sizeof((x)[0]))

#ifndef MIN
  #define MIN(a,b) (((a) < (b))? (a) : (b))
#endif

#ifndef MAX
  #define MAX(a,b) (((a) > (b))? (a) : (b))
#endif


//-----------------------------------------------------------------------------
// BASIC INTEGER TYPES

typedef  unsigned int  uint;

typedef  uint8_t   byte;

typedef  int8_t    int8;
typedef  int16_t   int16;
typedef  int32_t   int32;
typedef  int64_t   int64;

typedef  uint8_t   uint8;
typedef  uint16_t  uint16;
typedef  uint32_t  uint32;
typedef  uint64_t  uint64;


//-----------------------------------------------------------------------------
// BASIC FLOATING-POINT TYPES

typedef  float        float32;
typedef  double       float64;
typedef  long double  float80;
typedef  float80      real;

#define  FLOAT32_MANTISSA  23
#define  FLOAT64_MANTISSA  52
#define  FLOAT80_MANTISSA  68

#define  REAL_MANTISSA     FLOAT80_MANTISSA


//-----------------------------------------------------------------------------
// LINEAR INTERPOLATION

public_inline_function
real lerp(real t, real x0, real x1)
{
  //real x = ((1 - t) * x0) + (t * x1);
  real x = x0 + (t * (x1 - x0));
  return x;
}

public_inline_function
real unlerp(real x, real x0, real x1)
{
  real t = (x - x0) / (x1 - x0);
  return t;
}

public_inline_function
real swerp(real t, real x0, real x1)
{
  t = PI * t;        // Map [0,1] to [0,π] (i.e., 0° to 180°).
  t = cos(t);        // Map [0,π] to [1,-1].
  t = (1 - t) / 2;   // Map [1,-1] to [0,1].
  real x = lerp(t, x0, x1);
  return x;
}


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_function
void error_exit(char *message);


//-----------------------------------------------------------------------------
// PROJECT-LOCAL INCLUDE FILES

#import "MPReal.h"  // Wrapper for <mpfr.h>
#import "Vector3.h"
#import "Memory.h"
