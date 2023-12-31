/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Mandelbrot.h"


//-----------------------------------------------------------------------------
// NOTES ON ARBITRARY-PRECISION ARITHMETIC
//
// Comparison of runtimes:
//
// Area focusing on a baby Mandelbrot:
//
//     x_center:  -1.98999426899216
//     y_center:   0.00000000099994
//       x_size:   0.00000000000015
//     max_iter:   100000
//       pixels:   720 x 720
//     subpixel:   off
//
//     64-bit native hardware floating-point:     64.198 seconds
//     83-bit MPFR arbitrary precision:         5845.529 seconds (91x slower)
//
// Area with high noise and almost all exterior points (this is a very pretty
// picture, by the way):
//
//     x_center: -1.989994268992313
//     y_center:  0.000000000999976
//       x_size:  0.000000000000050
//     max_iter:  1000000
//       pixels:  see below
//     subpixel:  off
//
//     Size 32 x 32 pixels:
//     64-bit native hardware floating-point:      0.146 seconds
//     64-bit MPFR arbitrary precision:           10.343 seconds (71x slower)
//     128-bit MPFR arbitrary precision:          11.120 seconds (76x slower)
//     256-bit MPFR arbitrary precision:          12.962 seconds (89x slower)
//     512-bit MPFR arbitrary precision:          18.986 seconds (130x slower)
//     1024-bit MPFR arbitrary precision:         38.909 seconds (267x slower)
//     2048-bit MPFR arbitrary precision:         78.061 seconds (535x slower)
//     4096-bit MPFR arbitrary precision:        195.929 seconds (1342x slower)
//     8192-bit MPFR arbitrary precision:        530.959 seconds (3637x slower)
//
//     Size 256 x 256 pixels:
//     64-bit native hardware floating-point:      8.784 seconds
//     64-bit MPFR arbitrary precision:          614.287 seconds (70x slower)
//     128-bit MPFR arbitrary precision:         668.130 seconds (76x slower)
//     256-bit MPFR arbitrary precision:         771.506 seconds (88x slower)


//-----------------------------------------------------------------------------
// NOTES ON PERIODICITY CHECKING
//
// Calculations show that periodicity checking is only about 6% slower than
// brute-force iteration when calculating exterior points.
//
// For example, this deep, gnarly, non-interior region:
//
//    x_center: -1.98999426899005
//    y_center:  0.00000000100163
//      x_size:  0.00000000000150
//    max_iter:  100000 (chosen low enough to avoid discrimination between
//                       interior and exterior points)
//      pixels:  720x720
//    subpixel:  off
//
// has the following timings:
//
//    no periodicity checking:        39.018 seconds
//    exact periodicity checking:     41.296 seconds (5.84% slower)
//    inexact periodicity checking:   41.402 seconds (6.11% slower)
//
// For interior points, it is of course way faster -- anywhere from ten to
// thousands of times faster.


//-----------------------------------------------------------------------------
// ESCAPE/BAILOUT PARAMETERS

// An escape radius of 4 is ok, but 16 is much better for smooth coloring.
//private_static const real ESCAPE_RADIUS          = 16.0;
private_static const real ESCAPE_RADIUS_SQUARED  = 256.0;
private_static const real LOG2_ESCAPE_RADIUS     = 2.77258872223978;

// This adjusts the dwell value downward from the iteration value, because
// the dwell is supposed to be the number of iterations required to escape a
// radius of 2.  Note that this adjustment is not precise and is only a rough
// estimate based on (2^2)^2 = 16.
private_static const real ESCAPE_DWELL_OVERHEAD  = 2.0;

public_function
real mandelbrot_max_scalar_value_during_iteration()
{
  //return ESCAPE_RADIUS_SQUARED + 2.0;

  // FIXME: KLUDGE:  The value here is appropriate for the calculation
  // performed by the caller in Image.c, but it is not appropriate/correct
  // as a return value for this function.  To fix, rename the function to be
  // less specific.  The reason for using a low value here instead of the
  // actual value is because of the way floating-point arithmetic works.
  // Once escape is achieved from the radius-2 circle, the values are
  // guaranteed to escape to infinity, and any higher bits are unneeded in
  // precise calculations because there will be plenty already to the right
  // of the radix point.
  return 6.0;
}


//-----------------------------------------------------------------------------
// CREATE

public_constructor
Mandelbrot *mandelbrot_create(uint64 iter_max,
                              int mp_prec)
{
  assert(iter_max > 0);
  assert(mp_prec >= 16);

  Mandelbrot *this = mem_alloc(1, sizeof(*this));


  // Initialize configuration.

  this->conf.iter_max = iter_max;
  this->conf.mp_prec = mp_prec;


  // Initialize (zero out) statistics.

  this->stats.total_iter = 0;
  this->stats.total_probes = 0;
  this->stats.total_probes_uniterated = 0;

  this->stats.interior_iter = 0;
  this->stats.interior_probes = 0;
  this->stats.interior_probes_uniterated = 0;
  this->stats.interior_probes_aperiodic = 0;
  this->stats.interior_period_min = 0;
  this->stats.interior_period_max = 0;

  this->stats.exterior_iter = 0;
  this->stats.exterior_probes = 0;
  this->stats.exterior_probes_uniterated = 0;
  this->stats.exterior_dwell_min = NAN;
  this->stats.exterior_dwell_max = NAN;

  assert(ELEMENT_COUNT(this->stats.interior_probes_by_log2_iter) ==
         ELEMENT_COUNT(this->stats.exterior_probes_by_log2_iter));
  int k_max = ELEMENT_COUNT(this->stats.interior_probes_by_log2_iter) - 1;
  for (int k = 0; k <= k_max; k++)
  {
    this->stats.interior_probes_by_log2_iter[k] = 0;
    this->stats.exterior_probes_by_log2_iter[k] = 0;
  }


  return this;
}


//-----------------------------------------------------------------------------
// DESTROY

public_destructor
void mandelbrot_destroy(Mandelbrot **p_this)
{
  assert(p_this);
  Mandelbrot *this = *p_this;

  assert(this);

  mem_dealloc(p_this);
}


//-----------------------------------------------------------------------------
// COMPUTE DWELL FROM ITERATION COUNT AND FINAL Z VALUE

private_function
float64 dwell(uint64 iter, real z2)
{
  // FIXME: This calculation still isn't quite right.  For the image
  // (-.75,0),5 it produces fractional iterations in the range
  // -0.031193 to 0.971226, which has a width of 1.002188.
  // <http://en.wikipedia.org/wiki/Mandelbrot_set#Continuous_.28smooth.29_coloring>
  float64 fractional_iter = 1.5 - log2(log2(sqrt(z2)) / LOG2_ESCAPE_RADIUS);

  return iter + fractional_iter - ESCAPE_DWELL_OVERHEAD;
}


//-----------------------------------------------------------------------------
// CALCULATE NEXT ITERATION MARKER FOR PERIODICITY CHECKING
//
// NOTE to self:  Counting up the inner loop quadratically instead of
// exponentially (powers of 2) resulted in a 2x slowdown (for at least
// one tight deep area, anyway).  I had to try it because I saw someone
// else do it.  But in thinking more about it, it's no good because
// it means cycles can only be detected up to the square root of the
// maximum iteration count, as opposed to up to half of the maximum
// iteration count using the powers-of-2 method.  Exponential increase is
// the way to go.
//
// NOTE to self:  Using powers of 1.2 instead of 2 resulted in a 27%
// speedup.  I also tried a Fibonacci-sequence increase (which is powers
// of the golden ratio), but that only resulted in a 5% speedup over the
// powers of 2.  Then I tried 9/8, which was 33.4% better than 2.
// Then I tried 17/16, which was 35.9% better than 2.  Then I tried
// 101/100, which turned out to be 40.8% better than 2.  Wow.  Really
// slow exponential increases are the magic trick here.  I think I'll go
// with a hybrid approach for general use, because growth factors too
// close to 1 grow too slowly, which can be problematic for detecting
// high-period orbits, and because high growth factors are needed when
// the number of iterations is still low.
//
// In the table of intervals below, the left column (labeled "marker") contains
// values of the sequence produced by repeated calls to this routine, and the
// right column (labeled "cycle") contains the difference between successive
// values (in other words, the interval length), which is the maximum periodic
// orbit detectable by the interval.
//
// Marker Cycle    Marker Cycle    Marker Cycle    Marker Cycle    Marker Cycle
// ------ -----    ------ -----    ------ -----    ------ -----    ------ -----
//      1 1             2 2             4 4             8 8            16 8
//     24 12           36 18           54 27           81 20          101 25
//    126 31          157 39          196 49          245 61          306 38
//    344 43          387 48          435 54          489 61          550 68
//    618 77          695 86          781 97          878 109         987 123
//   1110 69         1179 73         1252 78         1330 83         1413 88
//   1501 93         1594 99         1693 105        1798 112        1910 119
//   2029 126        2155 134        2289 143        2432 152        2584 161
//   2745 171        2916 182        3098 193        3291 205        3496 218
//   3714 232        3946 246        4192 131        4323 135        4458 139
//   4597 143        4740 148        4888 152        5040 157        5197 162
//   5359 167        5526 172        5698 178        5876 183        6059 189
//   6248 195        6443 201        6644 207        6851 214        7065 220
//   7285 227        7512 234        7746 242        7988 249        8237 257
//   8494 265        8759 273        9032 282        9314 291        9605 300
//   9905 309       10214 319       10533 329       10862 339       11201 350
//  11551 360       11911 372       12283 383       12666 395       13061 408
//  13469 420       13889 434       14323 447       14770 461       15231 475
//  15706 490       16196 506       16702 521       17223 538       17761 555
//  18316 572       18888 590       19478 608       20086 627       20713 647
//  21360 667       22027 688       22715 709       23424 732       24156 754
//  24910 778       25688 802       26490 827       27317 853       28170 880
//  29050 907       29957 936       30893 965       31858 995       32853 1026
//  33879 1058      34937 1091      36028 1125      37153 1161      38314 1197
//  39511 1234      40745 1273      42018 1313      43331 1354      44685 1396
//  46081 1440      47521 1485      49006 1531      50537 1579      52116 1628
//  53744 1679      55423 1731      57154 1786      58940 1841      60781 1899
//  62680 1958      64638 2019      66657 2083      68740 2148      70888 2215
//  73103 2284      75387 2355      77742 2429      80171 2505      82676 2583
//  85259 2664      87923 2747      90670 2833      93503 2921      96424 3013
//  99437 3107     102544 3204     105748 3304     109052 3407     112459 3514
// 115973 3624     119597 3737     123334 3854     127188 3974     131162 4098
// 135260 4226     139486 4358     143844 4495     148339 4635     152974 4780
// 157754 4929     162683 5083     167766 5242     173008 5406     178414 5575
// 183989 5749     189738 5929     195667 6114     201781 6305     208086 6502
// 214588 6705     221293 6915     228208 7131     235339 7354     242693 7584
// 250277 7821     258098 8065     266163 8317     274480 8577     283057 8845
// 291902 9121     301023 9406     310429 9700     320129 10004    330133 10316
// 340449 10639    351088 10971    362059 11314    373373 11667    385040 12032
// 397072 12408    409480 12796    422276 13196    435472 13608    449080 14033
// 463113 14472    477585 14924    492509 15390    507899 15871    523770 16367
// 540137 16879    557016 17406    574422 17950    592372 18511    610883 19090
// 629973 19686    649659 20301    669960 20936    690896 21590    712486 22265
// 734751 22960    757711 23678    781389 24418    805807 25181    830988 25968
// 856956 26779    883735 27616    911351 28479    939830 29369    969199 30287

private_function
uint64 next_iteration_interval(uint64 iter, uint64 iter_max)
{
  uint64 grow = (iter <   16)?  1:  // Powers of  2/1  = 2
                (iter <   64)?  2:  // Powers of  3/2  = 1.5
                (iter <  256)?  4:  // Powers of  5/4  = 1.2
                (iter < 1024)?  8:  // Powers of  9/8  = 1.125
                (iter < 4096)? 16:  // Powers of 17/16 = 1.0625
                               32;  // Powers of 33/32 = 1.03125

  uint64 iter_bound = iter + MAX(iter/grow, 1);

  return MIN(iter_bound, iter_max);
}


//-----------------------------------------------------------------------------
// CALCULATE M-SET ITERATIONS FOR POINT

/*
  OBSOLETE -- This was the original periodicity-checking loop.
  NOTE to self:  Splitting the loop into two portions (thus avoiding having
  two separate conditionals on the iteration value inside the inner loop)
  results in a 33% speedup.  Way cool!

  real x = cx, y = cy, x_base = 1e99, y_base = 1e99, x2, y2;
  uint64 i, i_base = 0;
  for (i = 0; (i < i_max) && ((x2=x*x) + (y2=y*y) <= r2); i++)
  {
    if ((fabs((double)(x - x_base)) <= periodicity_epsilon) &&
        (fabs((double)(y - y_base)) <= periodicity_epsilon))
    {
      return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
    }
    else if (i == i_base)
    {
      x_base = x; y_base = y; i_base *= 2;
    }

    y = x * y; y += y + cy;
    x = x2 - y2 + cx;
  }

  if (i == i_max)
    return mandelbrot_result_interior_iterated_aperiodic(i_max);
  else
    return mandelbrot_result_exterior(i, dwell(i, x2 + y2));
*/


//-----------------------------------------------------------------------------
#if 0  // OBSOLETE
private_function
MandelbrotResult mandelbrot_compute_low_precision_periodicity_epsilon_old(
                   const real cx, const real cy, const real epsilon,
                   const uint64 i_max)
{
  assert(epsilon > 0);
  assert(i_max > 0);

  const real r2 = ESCAPE_RADIUS_SQUARED;

  real x = cx, y = cy, x_base = 1e99, y_base = 1e99;

  for (uint64 i_base = 0; i_base < i_max; )
  {
    uint64 i_bound = next_iteration_interval(i_base, i_max);

    for (uint64 i = i_base; i < i_bound; )
    {
      real x2 = x * x, y2 = y * y, z2 = x2 + y2;

      if (z2 > r2)
        return mandelbrot_result_exterior(i, dwell(i, z2));

      y = x * y; y += y + cy; x = x2 - y2 + cx; i++;

      if ((fabs((double)(x - x_base)) <= epsilon) &&
          (fabs((double)(y - y_base)) <= epsilon))
        return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
    }

    x_base = x; y_base = y; i_base = i_bound;
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}
#endif


//-----------------------------------------------------------------------------
// NOTE:  Unrolling the tight loop like this provides a 20% speedup over the
// unrolled version of this function.  Pretty cool.

private_function
MandelbrotResult mandelbrot_compute_low_precision_periodicity_epsilon(
                   const real cx, const real cy, const real epsilon,
                   const uint64 i_max)
{
  assert(epsilon > 0);
  assert(i_max > 0);

  // TODO:  Use square neighborhood when epsilon is less than or equal to 1,
  // and use round neighborhood when epsilon is greater than 1.  The reason
  // for this distinction is purely aesthetic.  Square is faster and should
  // always be used when the results of the check are invisible.  Round looks
  // better for large, sloppy epsilons used in low-quality quick renders.
  // The reason I'm not making a distiction here yet between square and round
  // is because it will require having two separate and nearly identical
  // iteration loops in order to maintain the efficiency of a single method.
  // So this is something for the future.
  #define EPSILON_SQUARE

  const real r2 = ESCAPE_RADIUS_SQUARED;
  #if !defined(EPSILON_SQUARE)
  const real epsilon2 = epsilon * epsilon;
  #endif

  real x = cx, y = cy, x_base = 1e99, y_base = 1e99;

  uint64 i = 0;
  uint64 i_base = i;
  uint64 i_bound = next_iteration_interval(i_base, i_max);
  while (i < i_max)
  {
    real save_x = x, save_y = y;
    real x2, y2, z2;
    bool cycle_detected = false;

    #define  ITERATE  \
      x2 = x * x; y2 = y * y; y = x * y; y += y + cy; x = x2 - y2 + cx;

    #if defined(EPSILON_SQUARE)  // Square neighborhood using L1 norm.
      #define  CHECK  \
        cycle_detected |= \
          (fabs((double)(x - x_base)) <= epsilon) && \
          (fabs((double)(y - y_base)) <= epsilon);
    #elif defined(EPSILON_ROUND)  // Round neighborhood using L2 norm.
      #define  CHECK  \
        cycle_detected |= \
          (fabs((double)(x - x_base)) <= epsilon) && \
          (fabs((double)(y - y_base)) <= epsilon) && \
          ((x-x_base)*(x-x_base) + (y-y_base)*(y-y_base) <= epsilon2)
    #elif defined(EPSILON_ROUND2)  // Round neighborhood using L2 norm.
      #define  CHECK  \
        cycle_detected |= \
          ((x-x_base)*(x-x_base) <= epsilon2) && \
          ((x-x_base)*(x-x_base) + (y-y_base)*(y-y_base) <= epsilon2)
    #endif
    // Sample timings from 
    //   cx=-.125 cy=.78 cd=0.33 i=1000000000 px=720 py=720
    //     ssimin=2 ssimax=4 ssemin=0 ssemax=4 sss=0.99
    // Square:  7484461749 iter / 54.970 sec = 136154186 iter/sec
    // Round1:  7794880719 iter / 59.116 sec = 131856884 iter/sec
    // Round2:  7794880719 iter / 70.539 sec = 110504502 iter/sec

    const int UNROLL = 8;
    ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK;
    ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK;
    i += UNROLL;

    #undef   ITERATE
    #undef   CHECK

    if ((z2 = x2 + y2) > 4)  // Use minimal escape radius here.
    {
      x = save_x; y = save_y; i -= UNROLL;
      while (i < i_max)
      {
        x2 = x * x; y2 = y * y;
        if ((z2 = x2 + y2) > r2)  // Use custom escape radius here.
          return mandelbrot_result_exterior_iterated(i+UNROLL, dwell(i, z2));
        y = x * y; y += y + cy; x = x2 - y2 + cx; i++;
      }
    }
    else if (cycle_detected)
    {
      return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
    }

    if (i >= i_bound)
    {
      x_base = x; y_base = y;
      i_base = i_bound; i_bound = next_iteration_interval(i_base, i_max);
    }
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}


//-----------------------------------------------------------------------------
#if 0  // OBSOLETE
private_function
MandelbrotResult mandelbrot_compute_low_precision_periodicity_exact_old(
                   const real cx, const real cy,
                   const uint64 i_max)
{
  assert(i_max > 0);

  const real r2 = ESCAPE_RADIUS_SQUARED;

  real x = cx, y = cy, x_base = 1e99, y_base = 1e99;

  for (uint64 i_base = 0; i_base < i_max; )
  {
    uint64 i_bound = next_iteration_interval(i_base, i_max);

    for (uint64 i = i_base; i < i_bound; )
    {
      real x2 = x * x, y2 = y * y, z2 = x2 + y2;

      if (z2 > r2)
        return mandelbrot_result_exterior(i, dwell(i, z2));

      y = x * y; y += y + cy; x = x2 - y2 + cx; i++;

      if ((x == x_base) && (y == y_base))
        return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
    }

    x_base = x; y_base = y; i_base = i_bound;
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}
#endif


//-----------------------------------------------------------------------------
// NOTE:  Unrolling the tight loop like this provides a 18% speedup over the
// unrolled version of this function.  Pretty cool.

private_function
MandelbrotResult mandelbrot_compute_low_precision_periodicity_exact(
                   const real cx, const real cy,
                   const uint64 i_max)
{
  assert(i_max > 0);

  const real r2 = ESCAPE_RADIUS_SQUARED;

  real x = cx, y = cy, x_base = 1e99, y_base = 1e99;

  uint64 i = 0;
  uint64 i_base = i;
  uint64 i_bound = next_iteration_interval(i_base, i_max);
  while (i < i_max)
  {
    real save_x = x, save_y = y;
    real x2, y2, z2;
    bool cycle_detected = false;

    #define  ITERATE  \
      x2 = x * x; y2 = y * y; y = x * y; y += y + cy; x = x2 - y2 + cx;

    #define  CHECK  \
      cycle_detected |= ((x == x_base) && (y == y_base));

    const int UNROLL = 8;
    ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK;
    ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK; ITERATE; CHECK;
    i += UNROLL;

    #undef   ITERATE
    #undef   CHECK

    if ((z2 = x2 + y2) > 4)  // Use minimal escape radius here.
    {
      x = save_x; y = save_y; i -= UNROLL;
      while (i < i_max)
      {
        x2 = x * x; y2 = y * y;
        if ((z2 = x2 + y2) > r2)  // Use custom escape radius here.
          return mandelbrot_result_exterior_iterated(i+UNROLL, dwell(i, z2));
        y = x * y; y += y + cy; x = x2 - y2 + cx; i++;
      }
    }
    else if (cycle_detected)
    {
      return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
    }

    if (i >= i_bound)
    {
      x_base = x; y_base = y;
      i_base = i_bound; i_bound = next_iteration_interval(i_base, i_max);
    }
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}


//-----------------------------------------------------------------------------
#if 0  // OBSOLETE
private_function
MandelbrotResult mandelbrot_compute_low_precision_no_periodicity_old(
                   const real cx, const real cy,
                   const uint64 i_max)
{
  assert(i_max > 0);

  const real r2 = ESCAPE_RADIUS_SQUARED;

  real x = cx, y = cy;

  for (uint64 i = 0; i < i_max; i++)
  {
    real x2 = x * x, y2 = y * y, z2 = x2 + y2;

    if (z2 > r2)
      return mandelbrot_result_exterior(i, dwell(i, z2));

    y = x * y; y += y + cy; x = x2 - y2 + cx;
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}
#endif


//-----------------------------------------------------------------------------
// NOTE:  Unrolling the tight loop like this provides an 18% speedup over the
// unrolled version.  Pretty cool.

private_function
MandelbrotResult mandelbrot_compute_low_precision_no_periodicity(
                   const real cx, const real cy,
                   const uint64 i_max)
{
  assert(i_max > 0);

  const real r2 = ESCAPE_RADIUS_SQUARED;

  real x = cx, y = cy;

  uint64 i = 0;
  while (i < i_max)
  {
    real save_x = x, save_y = y;
    real x2, y2, z2;

    #define  ITERATE  \
      x2 = x * x; y2 = y * y; y = x * y; y += y + cy; x = x2 - y2 + cx;

    const int UNROLL = 8;
    ITERATE; ITERATE; ITERATE; ITERATE; ITERATE; ITERATE; ITERATE; ITERATE;
    i += UNROLL;

    #undef  ITERATE

    if ((z2 = x2 + y2) > 4)  // Use minimal escape radius here.
    {
      x = save_x; y = save_y; i -= UNROLL;

      while (i < i_max)
      {
        x2 = x * x; y2 = y * y;
        if ((z2 = x2 + y2) > r2)  // Use custom escape radius here.
          return mandelbrot_result_exterior_iterated(i+UNROLL, dwell(i, z2));
        y = x * y; y += y + cy; x = x2 - y2 + cx; i++;
      }
    }
  }

  return mandelbrot_result_interior_iterated_aperiodic(i);
}


//-----------------------------------------------------------------------------
private_function
MandelbrotResult mandelbrot_compute_low_precision(
                   const real cx, const real cy,
                   const real periodicity_epsilon,
                   const uint64 iter_max)
{
  #if 0  // KLUDGE for testing -- makes a grid.
  if (((uint32)(int32)floor(cx) ^ (uint32)(int32)floor(cy)) & UINT32_C(1))
    return (MandelbrotResult) { .dwell = 0, .iter = 0, .period = 0 };
  else
    return mandelbrot_result_interior_uniterated();
  #endif

  // Early-out test for membership in main cardioid.
  if (cx >= -0.75)
  {
    const uint64 period = 1;
    real x = cx - 0.25;
    real y = x*x + cy*cy;
    x += y + y;
    if (x*x <= y)
      return mandelbrot_result_interior_uniterated_periodic(period);
  }

  // Early-out test for membership in largest disc, having center (-1,0) and
  // radius 1/4.
  else if (cx >= -1.25)
  {
    const uint64 period = 2;
    real x = cx + 1.0;
    if (x*x + cy*cy <= 0.0625)  // Radius 1/4 squared is 1/16.
      return mandelbrot_result_interior_uniterated_periodic(period);
  }

  // Early-out test for x-axis.  This is actually very, very important because
  // periodic orbits here are extremely rare and difficult to detect.  In fact,
  // in one case, before adding this test, I saw a region that chewed through
  // 7.89 trillion iterations (maximum iterations per probe was 1 billion).
  // After adding the test, the same image required only 18.3 million
  // iterations.  Thus, this simple little early-out test was quite worth it.
  // Generation of that image went from a whoppingly embarrassing 14.7 CPU
  // hours down to 1.9 CPU seconds!  This is probably the single most useful
  // optimization I have ever seen in my entire programming career.
  //
  // Note that the probe is recorded here as being aperiodic even though it
  // *might* be periodic.  This is not desirable, but the alternative (actually
  // determining whether the point is periodic) is just too expensive.
  else if (cx >= -2.0)
  {
    if (cy == 0)  // Note: It is already known that cx < -1.25.
      return mandelbrot_result_interior_uniterated_aperiodic();
  }

  // Handle other cases with either periodicity checking or standard counting.
  if (periodicity_epsilon > 0)
  {
    return mandelbrot_compute_low_precision_periodicity_epsilon(
                                    cx, cy, periodicity_epsilon, iter_max);
  }
  else if (periodicity_epsilon == 0)
  {
    return mandelbrot_compute_low_precision_periodicity_exact(cx, cy, iter_max);
  }
  else
  {
    return mandelbrot_compute_low_precision_no_periodicity(cx, cy, iter_max);
  }
}


//-----------------------------------------------------------------------------
private_function
MandelbrotResult mandelbrot_compute_high_precision(const mp_real cx,
                                                   const mp_real cy,
                                                   const mp_real epsilon,
                                                   const uint64 i_max)
{
  assert(i_max > 0);

  // Initialize constants and work variables.
  // FIXME:  These are never deallocated.
  static mp_real _n_0_75, _n_1_25, _n_2, _p_0_0625;
  static mp_real x, y, x2, y2, z2, r2, x_base, y_base, t;
  static bool initialized = false;
  if (!initialized)
  {
    mp_init2(_n_0_75, mp_get_prec(epsilon));
    mp_set_d(_n_0_75, -0.75);

    mp_init2(_n_1_25, mp_get_prec(epsilon));
    mp_set_d(_n_1_25, -1.25);

    mp_init2(_n_2, mp_get_prec(epsilon));
    mp_set_d(_n_2, -2.0);

    mp_init2(_p_0_0625, mp_get_prec(epsilon));
    mp_set_d(_p_0_0625, +0.0625);

    mp_init2(r2, mp_get_prec(epsilon));
    mp_set_d(r2, ESCAPE_RADIUS_SQUARED);

    mp_init2(x,      mp_get_prec(epsilon));
    mp_init2(y,      mp_get_prec(epsilon));
    mp_init2(x2,     mp_get_prec(epsilon));
    mp_init2(y2,     mp_get_prec(epsilon));
    mp_init2(z2,     mp_get_prec(epsilon));
    mp_init2(x_base, mp_get_prec(epsilon));
    mp_init2(y_base, mp_get_prec(epsilon));
    mp_init2(t,      mp_get_prec(epsilon));

    initialized = true;
  }

  // Early-out test for membership in main cardioid.
  if (mp_greaterequal_p(cx, _n_0_75))       // if (cx >= -0.75)
  {
    const uint64 period = 1;
    mp_sub_d(x, cx, 0.25);                  // x = cx - 0.25;
    mp_sqr(x2, x);                          //   x2 = x*x;
    mp_sqr(y2, cy);                         //   y2 = cy*cy;
    mp_add(y, x2, y2);                      // y = x*x + cy*cy;
    mp_add(x, x, y);                        // x += y;
    mp_add(x, x, y);                        // x += y;
    mp_sqr(x2, x);                          //   x2 = x*x;
    if (mp_lessequal_p(x2, y))              // if (x*x <= y)
      return mandelbrot_result_interior_uniterated_periodic(period);
  }

  // Early-out test for membership in largest disc.
  else if (mp_greaterequal_p(cx, _n_1_25))  // if (cx >= -1.25)
  {
    const uint64 period = 2;
    mp_add_d(x, cx, 1.0);                   // x = cx + 1.0;
    mp_sqr(x2, x);                          //   x2 = x*x;
    mp_sqr(y2, cy);                         //   y2 = cy*cy;
    mp_add(z2, x2, y2);                     //   z2 = x2 + y2;
    if (mp_lessequal_p(z2, _p_0_0625))      // if (x*x + cy*cy <= 0.0625)
      return mandelbrot_result_interior_uniterated_periodic(period);
  }

  // Early-out test for x-axis.
  else if (mp_greaterequal_p(cx, _n_2))     // if (cx >= -2.0)
  {
    // Important note:  The remainder of this test assumes that cx < -1.25,
    // having already been tested above.
    if (mp_zero_p(cy))
      return mandelbrot_result_interior_uniterated_aperiodic();
  }

  // Handle other cases with periodicity checking.

  mp_set(x, cx);                            // x = cx;
  mp_set(y, cy);                            // y = cy;
  mp_set_d(x_base, 1e99);                   // x_base = 1e99;
  mp_set_d(y_base, 1e99);                   // y_base = 1e99;
  // FIXME:  These ^^^^ are ugly and icky. Figure out some other way to set
  // a value of "infinity."

  for (uint64 i_base = 0; i_base < i_max; )
  {
    uint64 i_bound = next_iteration_interval(i_base, i_max);

    for (uint64 i = i_base; i < i_bound; )
    {
      mp_sqr(x2, x);                        // x2 = x * x;
      mp_sqr(y2, y);                        // y2 = y * y;
      mp_add(z2, x2, y2);                   // z2 = x2 + y2;

      if (mp_greater_p(z2, r2))             // if (z2 > r2)
      {
        real _z2 = mp_get_d(z2);          
        return mandelbrot_result_exterior_iterated(i, dwell(i, _z2));
      }

      mp_mul(y, x, y);                      // y = x * y;
      mp_add(y, y, y);                      // y += y;
      mp_add(y, y, cy);                     // y += cy;
      mp_sub(x, x2, y2);                    // x = x2 - y2;
      mp_add(x, x, cx);                     // x += cx;
      i++;

      mp_sub(t, x, x_base);                 // t = x - x_base;
      mp_abs(t, t);                         // t = fabs(x - x_base);
      if (mp_lessequal_p(t, epsilon))       // if (fabs(x-x_base) <= epsilon)
      {
        mp_sub(t, y, y_base);               // t = y - y_base;
        mp_abs(t, t);                       // t = fabs(y - y_base);
        if (mp_lessequal_p(t, epsilon))     // if (fabs(y-y_base) <= epsilon)
        {
          return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
        }
      }
    }

    mp_set(x_base, x);                      // x_base = x;
    mp_set(y_base, y);                      // y_base = y;
    i_base = i_bound;
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}


//-----------------------------------------------------------------------------
public_method
MandelbrotResult mandelbrot_compute(Mandelbrot *this,
                                    const mp_real cx, const mp_real cy,
                                    const mp_real periodicity_epsilon)
{
  assert(this);
  assert(this->conf.mp_prec >= mp_get_prec(periodicity_epsilon));
  assert(mp_sgn(periodicity_epsilon) >= 0);

  MandelbrotResult mr;


  // --- Compute sample.

  if (this->conf.mp_prec <= REAL_MANTISSA)
  {
    mr = mandelbrot_compute_low_precision(
           mp_get_d(cx),
           mp_get_d(cy),
           mp_get_d(periodicity_epsilon),
           this->conf.iter_max);
  }
  else
  {
    mr = mandelbrot_compute_high_precision(
           cx, cy,
           periodicity_epsilon,
           this->conf.iter_max);
  }


  // --- Accumulate statistics.

  this->stats.total_iter += mr.iter;
  this->stats.total_probes++;
  this->stats.total_probes_uniterated += (mr.iter == 0);

  if (mandelbrot_result_is_interior(mr))
  {
    this->stats.interior_iter += mr.iter;
    this->stats.interior_probes++;

    if (mandelbrot_result_is_interior_uniterated(mr))
    {
      this->stats.interior_probes_uniterated++;
    }
    else
    {
      assert(mandelbrot_result_is_interior_iterated(mr));

      if (mandelbrot_result_is_interior_aperiodic(mr))
      {
        this->stats.interior_probes_aperiodic++;
      }
      else
      {
        assert(mandelbrot_result_is_interior_periodic(mr));

        assert(mr.iter > 0);
        uint log2_iter = (uint)floor(log2(mr.iter));
        assert(log2_iter <
               ELEMENT_COUNT(this->stats.interior_probes_by_log2_iter));
        this->stats.interior_probes_by_log2_iter[log2_iter]++;
      }
    }

    if (mandelbrot_result_is_interior_periodic(mr))
       // Might be uniterated and that's okay!
    {
      if ((this->stats.interior_period_min == 0) ||
          (mr.period < this->stats.interior_period_min))
        this->stats.interior_period_min = mr.period;

      if ((this->stats.interior_period_max == 0) ||
          (mr.period > this->stats.interior_period_max))
        this->stats.interior_period_max = mr.period;
    }
  }
  else
  {
    this->stats.exterior_iter += mr.iter;
    this->stats.exterior_probes++;

    if (mandelbrot_result_is_exterior_uniterated(mr))
    {
      this->stats.exterior_probes_uniterated++;
    }
    else
    {
      assert(mr.iter > 0);
      uint log2_iter = (uint)floor(log2(mr.iter));
      assert(log2_iter <
             ELEMENT_COUNT(this->stats.exterior_probes_by_log2_iter));
      this->stats.exterior_probes_by_log2_iter[log2_iter]++;
    }

    if (isnan(this->stats.exterior_dwell_min) ||
        (mr.dwell < this->stats.exterior_dwell_min))
      this->stats.exterior_dwell_min = mr.dwell;

    if (isnan(this->stats.exterior_dwell_max) ||
        (mr.dwell > this->stats.exterior_dwell_max))
      this->stats.exterior_dwell_max = mr.dwell;
  }


  return mr;
}


//-----------------------------------------------------------------------------
