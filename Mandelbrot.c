/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Mandelbrot.h"


//-----------------------------------------------------------------------------
// BAILOUT PARAMETERS

// A bailout radius of 4 is ok, but 16 is much better for smooth coloring.
private_static const real BAILOUT_RADIUS          = 16.0;
private_static const real BAILOUT_RADIUS_SQUARED  = 256.0;
private_static const real LOG_BAILOUT_RADIUS      = 2.77258872223978;

// This adjusts the dwell value downward from the iteration value, because
// the dwell is supposed to be the number of iterations required to escape a
// radius of 2.  Note that this adjustment is not precise and is only a rough
// estimate based on (2^2)^2 = 16.
private_static const real BAILOUT_DWELL_OVERHEAD  = 2.0;


//-----------------------------------------------------------------------------
// COMPUTE DWELL FROM ITERATION COUNT AND FINAL Z VALUE

private_function
float64 dwell(uint64 iter, real z2)
{
  // FIXME: This calculation still isn't quite right.  For the image
  // (-.75,0),5 it produces fractional iterations in the range
  // -0.031193 to 0.971226, which has a width of 1.002188.
  // <http://en.wikipedia.org/wiki/Mandelbrot_set#Continuous_.28smooth.29_coloring>
  float64 fractional_iter = 1.5 - log2(log2(sqrt(z2)) / LOG_BAILOUT_RADIUS);

  return iter + fractional_iter - BAILOUT_DWELL_OVERHEAD;
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
  for (i = 0; (i < i_max) && ((x2=x*x) + (y2=y*y) < r2); i++)
  {
    if ((fabs(x - x_base) <= periodicity_epsilon) &&
        (fabs(y - y_base) <= periodicity_epsilon))
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
private_function
MandelbrotResult mandelbrot_compute_with_periodicity_epsilon(
                   const real cx, const real cy,
                   const real epsilon,
                   const uint64 i_max)
{
  assert(epsilon > 0);
  assert(i_max > 0);

  const real r2 = BAILOUT_RADIUS_SQUARED;

  real x = cx, y = cy, x_base = 1e99, y_base = 1e99;

  for (uint64 i_base = 0; i_base < i_max; )
  {
    uint64 i_bound = next_iteration_interval(i_base, i_max);

    for (uint64 i = i_base; i < i_bound; )
    {
      real x2 = x * x, y2 = y * y, z2 = x2 + y2;

      if (z2 >= r2)
        return mandelbrot_result_exterior(i, dwell(i, z2));

      y = x * y; y += y + cy; x = x2 - y2 + cx; i++;

      if ((fabs(x - x_base) <= epsilon) && (fabs(y - y_base) <= epsilon))
        return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
    }

    x_base = x; y_base = y; i_base = i_bound;
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}


//-----------------------------------------------------------------------------
private_function
MandelbrotResult mandelbrot_compute_with_periodicity_exact(const real cx,
                                                           const real cy,
                                                           const uint64 i_max)
{
  assert(i_max > 0);

  const real r2 = BAILOUT_RADIUS_SQUARED;

  real x = cx, y = cy, x_base = 1e99, y_base = 1e99;

  for (uint64 i_base = 0; i_base < i_max; )
  {
    uint64 i_bound = next_iteration_interval(i_base, i_max);

    for (uint64 i = i_base; i < i_bound; )
    {
      real x2 = x * x, y2 = y * y, z2 = x2 + y2;

      if (z2 >= r2)
        return mandelbrot_result_exterior(i, dwell(i, z2));

      y = x * y; y += y + cy; x = x2 - y2 + cx; i++;

      if ((x == x_base) && (y == y_base))
        return mandelbrot_result_interior_iterated_periodic(i, i - i_base);
    }

    x_base = x; y_base = y; i_base = i_bound;
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}


//-----------------------------------------------------------------------------
private_function
MandelbrotResult mandelbrot_compute_with_no_periodicity(const real cx,
                                                        const real cy,
                                                        const uint64 i_max)
{
  assert(i_max > 0);

  const real r2 = BAILOUT_RADIUS_SQUARED;

  real x = cx, y = cy;

  for (uint64 i = 0; i < i_max; i++)
  {
    real x2 = x * x, y2 = y * y, z2 = x2 + y2;

    if (z2 >= r2)
      return mandelbrot_result_exterior(i, dwell(i, z2));

    y = x * y; y += y + cy; x = x2 - y2 + cx;
  }

  return mandelbrot_result_interior_iterated_aperiodic(i_max);
}


//-----------------------------------------------------------------------------
public_function
MandelbrotResult mandelbrot_compute(const real cx, const real cy,
                                    const real periodicity_epsilon,
                                    const uint64 iter_max)
{
  assert(iter_max > 0);

  // Early-out test for membership in largest disc.
  if (cx <= -0.75)
  {
    real x = cx + 1.0;
    if (x*x + cy*cy < 0.0625)
      return mandelbrot_result_interior_uniterated();
  }

  // Early-out test for membership in main cardioid.
  else
  {
    real x = cx - 0.25;
    real y = x*x + cy*cy;
    x += y + y;
    if (x*x <= y)
      return mandelbrot_result_interior_uniterated();
  }

  // Handle other cases with either periodicity checking or standard counting.
  // Calculations show that periodicity checking is only about 15% slower than
  // brute-force iteration when calculating exterior points.  For interior
  // points, it is of course way faster -- anywhere from 10 to thousands of
  // times faster.
  if (periodicity_epsilon > 0)
  {
    return mandelbrot_compute_with_periodicity_epsilon(
                                    cx, cy, periodicity_epsilon, iter_max);
  }
  else if (periodicity_epsilon == 0)
  {
    return mandelbrot_compute_with_periodicity_exact(cx, cy, iter_max);
  }
  else
  {
    return mandelbrot_compute_with_no_periodicity(cx, cy, iter_max);
  }
}


//-----------------------------------------------------------------------------
