/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "DeepReal.h"

//#include <gmp.h>
#include <mpfr.h>


//-----------------------------------------------------------------------------
// USEFUL CONSTANTS

// 2^32
#define  DEEP_RADIX  4294967296.0

#if 0  // OBSOLETE
#define  SIGN_BIT               (UINT32_C(1) << 31)
#define  SIGN_BIT_IS_SET(this)  ((this->v[0] &   SIGN_BIT) != 0)
#define  SET_SIGN_BIT(this)       this->v[0] |=  SIGN_BIT
#define  CLEAR_SIGN_BIT(this)     this->v[0] &= ~SIGN_BIT
#endif


//-----------------------------------------------------------------------------
// ALLOCATE

public_method
DeepReal *dr_alloc(int bits_deep)
{
  assert(bits_deep >= 32);

  int bits_total = 32 + bits_deep;
  int element_count = (bits_total + 31) / 32;

  DeepReal *this = mem_alloc(1, sizeof(DeepReal) +
                                (element_count * sizeof(this->v[0])));
  this->k = element_count;
  for (int i = 0; i < this->k; i++)
    this->v[i] = 0;

  return this;
}


//-----------------------------------------------------------------------------
// ALLOCATE WITH INTEGER VALUE

public_method
DeepReal *dr_alloc_int32(int bits_deep, int32 n)
{
  DeepReal *this = dr_alloc(bits_deep);
  dr_set_int32(this, n);
  return this;
}


//-----------------------------------------------------------------------------
// ALLOCATE WITH FLOATING-POINT VALUE

public_method
DeepReal *dr_alloc_float64(int bits_deep, float64 x)
{
  DeepReal *this = dr_alloc(bits_deep);
  dr_set_float64(this, x);
  return this;
}


//-----------------------------------------------------------------------------
// DEALLOCATE

public_method
void dr_dealloc(DeepReal **p_this)
{
  assert(p_this);
  DeepReal *this = *p_this;

  assert(this);

  mem_dealloc(p_this);
}


//-----------------------------------------------------------------------------
public_method
DeepReal *dr_clone(DeepReal *that)
{
  assert(that);

  return (DeepReal*)mem_alloc_clone(that);
}


//-----------------------------------------------------------------------------
public_method
void dr_print_debug(DeepReal *this)
{
  assert(this);

  printf("{s=%d, k=%d, v=[", (int)this->s, this->k);
  for (int i = 0; i < this->k; i++)
  {
    printf("%08X%s", this->v[i], (i < this->k - 1)? " ":"");
  }
  printf("]}");
}


//-----------------------------------------------------------------------------
public_method
void dr_print(DeepReal *this)
{
  assert(this);

  DeepReal *temp = dr_clone(this);

  printf("%s%u.", temp->s? "-":"+", temp->v[0]);
  temp->v[0] = 0;

  int base10_digits_deep = floor((double)(32 * (temp->k - 1)) / log2(10));
  for (int i = 0; i < base10_digits_deep; i++)
  {
    dr_mul_int32(temp, 10);
    int digit = temp->v[0];
    assert((digit >= 0) && (digit <= 9));
    printf("%d", digit);
    temp->v[0] = 0;
  }

  dr_dealloc(&temp);
}


//-----------------------------------------------------------------------------
public_method
void dr_assign(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  this->s = that->s;
  this->k = that->k;
  mem_copy(this->v, that->v, this->k, sizeof(this->v[0]));
}


//-----------------------------------------------------------------------------
public_method
void dr_negate(DeepReal *this)
{
  assert(this);

  this->s = !this->s;
}


//-----------------------------------------------------------------------------
public_method
float64 dr_get_float64(DeepReal *this)
{
  assert(this);

  float64 x = 0;

  for (int i = this->k - 1; i >= 0; i--)
  {
    x /= DEEP_RADIX;
    x += this->v[i];
  }

  if (this->s)
    x = -x;

  return x;
}


//-----------------------------------------------------------------------------
public_method
void dr_set_float64(DeepReal *this, float64 x)
{
  assert(this);
  assert((x > -(DEEP_RADIX/2)) && (x < +(DEEP_RADIX/2)));

  if (x < 0)  { this->s = true;  x = -x; }
  else        { this->s = false;         }

  for (int i = 0; i < this->k; i++)
  {
    this->v[i] = (uint32)floor(x);
    x -= (float64)this->v[i];
    x *= DEEP_RADIX;
  }
}


//-----------------------------------------------------------------------------
public_method
int32 dr_get_int32(DeepReal *this)
{
  assert(this);

  if (this->s)
    return -(int32)this->v[0];
  else
    return +(int32)this->v[0];
}


//-----------------------------------------------------------------------------
public_method
void dr_set_int32(DeepReal *this, int32 n)
{
  assert(this);

  if (n < 0)  { this->s = true;  n = -n; }
  else        { this->s = false;         }

  this->v[0] = n;
  for (int i = 1; i < this->k; i++)
    this->v[i] = 0;
}


//-----------------------------------------------------------------------------
public_method
void dr_mul_two(DeepReal *this)
{
  assert(this);

  uint32 ci = 0;  // Carry-in bit from lower segment.
  for (int i = this->k - 1; i >= 0; i--)
  {
    uint32 co = this->v[i] >> 31;  // Carry-out bit to higher segment.
    this->v[i] = (this->v[i] << 1) | ci;
    ci = co;
  }
}


//-----------------------------------------------------------------------------
public_method
void dr_mul_int32(DeepReal *this, int32 n)
{
  assert(this);

  if (n < 0)  { this->s = !this->s; n = -n; }

  uint32 c = 0;  // Carry segment.
  for (int i = this->k - 1; i >= 0; i--)
  {
    uint64 p = (uint64)c + ((uint64)(uint32)n * (uint64)this->v[i]); // Product.
    this->v[i] = (uint32)p;
    c = (uint32)(p >> 32);
  }
}


//-----------------------------------------------------------------------------
public_method
void dr_div_two(DeepReal *this)
{
  assert(this);

  uint32 v_bor = 0;  // Bitwise-or of all segments (will be computed).
  uint32 ci = 0;  // Carry-in bit from higher segment.
  for (int i = 0; i < this->k; i++)
  {
    uint32 co = this->v[i] & 1;  // Carry-out bit to lower segment.
    this->v[i] = (this->v[i] >> 1) | (ci << 31);
    v_bor |= this->v[i];
    ci = co;
  }

  // Round the quotient.  Increment last bit (and propagate carry) if remainder
  // of last division is more than half the divisor.  This looks like another
  // loop, but in practice it will only go one time around in the loop in the
  // vast majority of cases.  NOTE:  This step is *skipped* if the final result
  // is zero, to avoid 0.0000000000000000...1 mapping back to itself.
  if (ci)
  {
    unless (v_bor == 0)
    {
      for (int i = this->k - 1; i >= 0; i--)
        if ((this->v[i] += 1) != 0)
          break;
    }
  }
}


//-----------------------------------------------------------------------------
public_method
void dr_div_int32(DeepReal *this, int32 n)
{
  assert(this);
  assert(n != 0);

  if (n < 0)  { this->s = !this->s; n = -n; }

  uint32 v_bor = 0;  // Bitwise-or of all segments (will be computed).
  uint32 r = 0;  // Remainder.
  for (int i = 0; i < this->k; i++)
  {
    uint64 d = ((uint64)r << 32) + (uint64)this->v[i];

    this->v[i] = (uint32)(d / (uint64)(uint32)n);
    r          = (uint32)(d % (uint64)(uint32)n);

    v_bor |= this->v[i];
  }

  // Round the quotient.  Increment last bit (and propagate carry) if remainder
  // of last division is more than half the divisor.  This looks like another
  // loop, but in practice it will only go one time around in the loop in the
  // vast majority of cases.  NOTE:  This step is *skipped* in the very special
  // case where the final result is zero and the divisor is two, in order to
  // avoid 0.0000000000000000...1 mapping back to itself.
  if (2 * r >= (uint32)n)
  {
    unless ((v_bor == 0) && (n == 2))
    {
      for (int i = this->k - 1; i >= 0; i--)
        if ((this->v[i] += 1) != 0)
          break;
    }
  }
}


//-----------------------------------------------------------------------------
private_method
void dr_normalize_sign_internal(DeepReal *this)
{
  assert(this);

  if (this->v[0] & (UINT32_C(1) << 31))
  {
    uint32 c = 1;
    for (int i = this->k - 1; i >= 0; i--)
    {
      uint64 s = (uint64)(~(this->v[i])) + (uint64)c;
      this->v[i] = (uint32)s;
      c = (uint32)(s >> 32);
    }

    this->s = !this->s;
  }
}


//-----------------------------------------------------------------------------
private_method
void dr_add_unsigned_internal(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  uint32 c = 0;
  for (int i = this->k - 1; i >= 0; i--)
  {
    uint64 s = (uint64)this->v[i] + (uint64)that->v[i] + (uint64)c;
    this->v[i] = (uint32)s;
    c = (uint32)(s >> 32);
  }
}


//-----------------------------------------------------------------------------
private_method
void dr_sub_unsigned_internal(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  uint32 b = 0;
  for (int i = this->k - 1; i >= 0; i--)
  {
    uint64 d = (uint64)this->v[i] - (uint64)that->v[i] - (uint64)b;
    this->v[i] = (uint32)d;
    b = (uint32)(d >> 32);
  }
}


//-----------------------------------------------------------------------------
public_method
void dr_add(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (!this->s && !that->s)
  {
    // Positive plus Positive
    dr_add_unsigned_internal(this, that);
  }
  else if (!this->s && that->s)
  {
    // Positive plus Negative
    dr_sub_unsigned_internal(this, that);
    dr_normalize_sign_internal(this);
  }
  else if (this->s && !that->s)
  {
    // Negative plus Positive
    dr_sub_unsigned_internal(this, that);
    dr_normalize_sign_internal(this);
  }
  else if (this->s && that->s)
  {
    // Negative plus Negative
    dr_add_unsigned_internal(this, that);
  }
}


//-----------------------------------------------------------------------------
public_method
void dr_sub(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (!this->s && !that->s)
  {
    // Positive minus Positive
    dr_sub_unsigned_internal(this, that);
    dr_normalize_sign_internal(this);
  }
  else if (!this->s && that->s)
  {
    // Positive minus Negative
    dr_add_unsigned_internal(this, that);
  }
  else if (this->s && !that->s)
  {
    // Negative minus Positive
    dr_add_unsigned_internal(this, that);
  }
  else if (this->s && that->s)
  {
    // Negative minus Negative
    dr_sub_unsigned_internal(this, that);
    dr_normalize_sign_internal(this);
  }
}


//-----------------------------------------------------------------------------
public_method
bool dr_eq(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (this->s != that->s)
    return false;
  for (int i = 0; i < this->k; i++)
    if (this->v[i] != that->v[i])
      return false;
  return true;
}


//-----------------------------------------------------------------------------
public_method
bool dr_ne(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (this->s != that->s)
    return true;
  for (int i = 0; i < this->k; i++)
    if (this->v[i] != that->v[i])
      return true;
  return false;
}


//-----------------------------------------------------------------------------
private_method
bool dr_lt_unsigned_internal(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  for (int i = 0; i < this->k; i++)
  {
         if (this->v[i] < that->v[i]) return true;
    else if (this->v[i] > that->v[i]) return false;
    else                              continue;
  }
  return false;
}


//-----------------------------------------------------------------------------
private_method
bool dr_lte_unsigned_internal(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  for (int i = 0; i < this->k; i++)
  {
         if (this->v[i] < that->v[i]) return true;
    else if (this->v[i] > that->v[i]) return false;
    else                              continue;
  }
  return true;
}


//-----------------------------------------------------------------------------
private_method
bool dr_gt_unsigned_internal(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  for (int i = 0; i < this->k; i++)
  {
         if (this->v[i] > that->v[i]) return true;
    else if (this->v[i] < that->v[i]) return false;
    else                              continue;
  }
  return false;
}


//-----------------------------------------------------------------------------
private_method
bool dr_gte_unsigned_internal(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  for (int i = 0; i < this->k; i++)
  {
         if (this->v[i] > that->v[i]) return true;
    else if (this->v[i] < that->v[i]) return false;
    else                              continue;
  }
  return true;
}


//-----------------------------------------------------------------------------
public_method
bool dr_lt(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (this->s && !that->s)
  {
    // Negative compared to Positive
    return true;
  }
  else if (this->s && that->s)
  {
    // Negative compared to Negative
    return dr_gt_unsigned_internal(this, that);
  }
  else if (!this->s && !that->s)
  {
    // Positive compared to Positive
    return dr_lt_unsigned_internal(this, that);
  }
  else if (!this->s && that->s)
  {
    // Positive compared to Negative
    return false;
  }
  assert(false);
  return false;
}


//-----------------------------------------------------------------------------
public_method
bool dr_lte(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (this->s && !that->s)
  {
    // Negative compared to Positive
    return true;
  }
  else if (this->s && that->s)
  {
    // Negative compared to Negative
    return dr_gte_unsigned_internal(this, that);
  }
  else if (!this->s && !that->s)
  {
    // Positive compared to Positive
    return dr_lte_unsigned_internal(this, that);
  }
  else if (!this->s && that->s)
  {
    // Positive compared to Negative
    return false;
  }
  assert(false);
  return false;
}


//-----------------------------------------------------------------------------
public_method
bool dr_gt(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (this->s && !that->s)
  {
    // Negative compared to Positive
    return false;
  }
  else if (this->s && that->s)
  {
    // Negative compared to Negative
    return dr_lt_unsigned_internal(this, that);
  }
  else if (!this->s && !that->s)
  {
    // Positive compared to Positive
    return dr_gt_unsigned_internal(this, that);
  }
  else if (!this->s && that->s)
  {
    // Positive compared to Negative
    return true;
  }
  assert(false);
  return false;
}

//-----------------------------------------------------------------------------
public_method
bool dr_gte(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that); assert(this->k == that->k);

  if (this->s && !that->s)
  {
    // Negative compared to Positive
    return false;
  }
  else if (this->s && that->s)
  {
    // Negative compared to Negative
    return dr_lte_unsigned_internal(this, that);
  }
  else if (!this->s && !that->s)
  {
    // Positive compared to Positive
    return dr_gte_unsigned_internal(this, that);
  }
  else if (!this->s && that->s)
  {
    // Positive compared to Negative
    return true;
  }
  assert(false);
  return false;
}


//-----------------------------------------------------------------------------
public_method
void dr_mul(DeepReal *this, DeepReal *that)
{
  assert(this); assert(that);
}


//-----------------------------------------------------------------------------
public_method
void dr_sqr(DeepReal *this)
{
  assert(this);
}


//-----------------------------------------------------------------------------
public_function
void dr_test(void)
{
#if 0
  int bits_deep = 32*10;
  DeepReal *x_min = dr_alloc(bits_deep);
  DeepReal *x_max = dr_alloc(bits_deep);
  DeepReal *dx    = dr_alloc(bits_deep);
  DeepReal *x     = dr_alloc(bits_deep);
  DeepReal *zero  = dr_alloc(bits_deep);

  dr_set_float64(x_min, -2.00);
  dr_set_float64(x_max,  0.50);
  dr_assign(dx, x_max); dr_sub(dx, x_min);
  dr_set_int32(zero, 0);

  int dj = 49;
  for (int j = 0; j < dj; j++)
  {
    dr_assign(x, dx);
    dr_mul_int32(x, j);
    dr_div_int32(x, dj - 1);
    dr_add(x, x_min);

    printf("%4d.  ", j);
    dr_print(x);
    printf("\n");
  }

  printf("x_min = "); dr_print(x_min); printf("\n");
  printf("x_max = "); dr_print(x_max); printf("\n");
  printf("dx    = "); dr_print(dx);    printf("\n");

  printf("\n\n");
  dr_set_float64(x, 2.5);
  dr_mul_int32(x, 3);
  dr_div_int32(x, 48);
  dr_set_int32(x, 1);
  while (dr_ne(x, zero))
  {
    dr_mul_int32(x, 2);
    //dr_print(x); printf("\n");
    dr_div_int32(x, 5);
    //dr_print(x); printf("\n");
  }
  dr_print(x); printf("\n");
#endif


  // Hmmm... Well... For this test below, which computes e to high precision,
  // my code takes 1.798s, and MPFR takes only 0.482s -- so MPFR is almost 4
  // times faster than my implementation here.  This isn't too surprising,
  // since it is a very mature library, and it uses true 64-bit operations
  // beneath the hood, and most importantly it uses hand-coded/hand-tuned
  // assembly language.  So I think I will abandon my approach and just use
  // MPFR instead.

  int bits_deep = 65536; //32*10;
#if 1
  DeepReal *a = dr_alloc(bits_deep);
  DeepReal *x = dr_alloc(bits_deep);
  dr_set_int32(a, 0);
  dr_set_int32(x, 1);
  for (int i = 1; i <= bits_deep; i++)
  {
    dr_add(a, x);
    dr_div_int32(x, i);
  }
  dr_print(a); printf("\n");
#else
  mpfr_t a, x;
  mpfr_init2(a, bits_deep);
  mpfr_init2(x, bits_deep);
  mpfr_set_d(a, 0, MPFR_RNDN);
  mpfr_set_d(x, 1, MPFR_RNDN);
  for (int i = 1; i <= bits_deep; i++)
  {
    mpfr_add(a, a, x, MPFR_RNDN);
    mpfr_div_ui(x, x, i, MPFR_RNDN);
  }
  mpfr_out_str(stdout, 10, 0, a, MPFR_RNDN);
  printf("\n");
  mpfr_clear(x);
  mpfr_clear(a);
#endif
}


//-----------------------------------------------------------------------------
