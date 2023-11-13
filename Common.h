/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

// Standard library
#import <stdbool.h>
#import <stdint.h>
#import <stdlib.h>
#import <stdio.h>
#import <unistd.h>
#import <memory.h>
#import <malloc/malloc.h>
#import <string.h>
#import <ctype.h>
#import <math.h>
#import <complex.h>
#import <limits.h>
#import <assert.h>

// Additional libraries
#import <mpfr.h>


//------------------------------------------------------------------------------
// The overloading of the word "static" in C is annoying.  Sometimes it actually
// means static, but many times it means private, and sometimes it is even
// needed for public functions (namely, when declaring inline functions in a
// header file).  For easier reading, the following are used in this program
// instead.
//
#define  private_static            static
#define  public_static

#define  private_const             static const
#define  public_const              static const

#define  private_function          static
#define  public_function
#define  extern_public_function    extern

#define  private_method            static
#define  public_method  
#define  extern_public_method      extern

#define  private_inline_function   static inline
#define  public_inline_function    static inline

#define  private_inline_method     static inline
#define  public_inline_method      static inline


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

typedef  double  real;

typedef  float   float32;
typedef  double  float64;


//-----------------------------------------------------------------------------
// LINEAR INTERPOLATION

public_inline_function
double lerp(double t, double x0, double x1)
{
  //double x = ((1 - t) * x0) + (t * x1);
  double x = x0 + (t * (x1 - x0));
  return x;
}

public_inline_function
double unlerp(double x, double x0, double x1)
{
  double t = (x - x0) / (x1 - x0);
  return t;
}

public_inline_function
double swerp(double t, double x0, double x1)
{
  t = PI * t;        // Map [0,1] to [0,π] (i.e., 0° to 180°).
  t = cos(t);        // Map [0,π] to [1,-1].
  t = (1 - t) / 2;   // Map [1,-1] to [0,1].
  double x = lerp(t, x0, x1);
  return x;
}


//-----------------------------------------------------------------------------
// FUNCTION PROTOTYPES

extern_public_function
void error_exit(char *message);


//-----------------------------------------------------------------------------
// PROJECT-LOCAL INCLUDE FILES

#import "Memory.h"
