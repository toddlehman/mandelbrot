/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//------------------------------------------------------------------------------
// Include MPFR multiprecision floating-point library and modify its interface
// to be a bit cleaner.

#if USE_MPFR  // (Manifest constant defined in Makefile)

#import <mpfr.h>


// I hate the "_t" suffix for "type".  Hate it, hate it, hate it.
typedef  mpfr_t  mp_real;


// This program always uses MPFR_RNDN.

#define  MP_ROUND  MPFR_RNDN


// Shorten all the function names.

#define  mp_init                       mpfr_init
#define  mp_init2                      mpfr_init2
#define  mp_clear                      mpfr_clear
#define  mp_clears                     mpfr_clears

#define  mp_get_default_prec           mpfr_get_default_prec
#define  mp_set_default_prec           mpfr_set_default_prec
#define  mp_get_prec                   mpfr_get_prec
#define  mp_set_prec                   mpfr_set_prec

#define  mp_set(a,b)                   mpfr_set(a,b,MP_ROUND)
#define  mp_set_ui(a,b)                mpfr_set_ui(a,b,MP_ROUND)
#define  mp_set_si(a,b)                mpfr_set_si(a,b,MP_ROUND)
#define  mp_set_uj(a,b)                mpfr_set_uj(a,b,MP_ROUND)
#define  mp_set_sj(a,b)                mpfr_set_sj(a,b,MP_ROUND)
#define  mp_set_flt(a,b)               mpfr_set_flt(a,b,MP_ROUND)
#define  mp_set_d(a,b)                 mpfr_set_d(a,b,MP_ROUND)
#define  mp_set_ld(a,b)                mpfr_set_ld(a,b,MP_ROUND)
#define  mp_set_decimal64(a,b)         mpfr_set_decimal64(a,b,MP_ROUND)
#define  mp_set_z(a,b)                 mpfr_set_z(a,b,MP_ROUND)
#define  mp_set_q(a,b)                 mpfr_set_q(a,b,MP_ROUND)
#define  mp_set_f(a,b)                 mpfr_set_f(a,b,MP_ROUND)

#define  mp_set_ui_2exp(a,b,c)         mpfr_set_ui_2exp(a,b,c,MP_ROUND)
#define  mp_set_si_2exp(a,b,c)         mpfr_set_si_2exp(a,b,c,MP_ROUND)
#define  mp_set_uj_2exp(a,b,c)         mpfr_set_uj_2exp(a,b,c,MP_ROUND)
#define  mp_set_sj_2exp(a,b,c)         mpfr_set_sj_2exp(a,b,c,MP_ROUND)
#define  mp_set_z_2exp(a,b,c)          mpfr_set_z_2exp(a,b,c,MP_ROUND)

#define  mp_set_str(a,b,c)             mpfr_set_str(a,b,c,MP_ROUND)

#define  mp_strtofr(a,b,c,d)           mpfr_strtofr(a,b,c,d,MP_ROUND)

#define  mp_set_nan                    mpfr_set_nan
#define  mp_set_inf                    mpfr_set_inf
#define  mp_set_zero                   mpfr_set_zero

#define  mp_swap                       mpfr_swap

#define  mp_init_set                   mpfr_init_set
#define  mp_init_set_ui                mpfr_init_set_ui
#define  mp_init_set_si                mpfr_init_set_si
#define  mp_init_set_d                 mpfr_init_set_d
#define  mp_init_set_ld                mpfr_init_set_ld
#define  mp_init_set_z                 mpfr_init_set_z
#define  mp_init_set_q                 mpfr_init_set_q
#define  mp_init_set_f                 mpfr_init_set_f
#define  mp_init_set_str               mpfr_init_set_str

#define  mp_get_flt(a)                 mpfr_get_flt(a,MP_ROUND)
#define  mp_get_d(a)                   mpfr_get_d(a,MP_ROUND)
#define  mp_get_ld(a)                  mpfr_get_ld(a,MP_ROUND)
#define  mp_get_decimal64(a)           mpfr_get_decimal64(a,MP_ROUND)

#define  mp_get_si(a)                  mpfr_get_si(a,MP_ROUND)
#define  mp_get_ui(a)                  mpfr_get_ui(a,MP_ROUND)
#define  mp_get_sj(a)                  mpfr_get_sj(a,MP_ROUND)
#define  mp_get_uj(a)                  mpfr_get_uj(a,MP_ROUND)

#define  mp_get_d_2exp(a,b)            mpfr_get_d_2exp(a,b,MP_ROUND)
#define  mp_get_ld_2exp(a,b)           mpfr_get_ld_2exp(a,b,MP_ROUND)

#define  mp_frexp(a,b,c)               mpfr_frexp(a,b,c,MP_ROUND)

#define  mp_get_z_2exp                 mpfr_get_z_2exp

#define  mp_get_z(a,b,c)               mpfr_get_z(a,b,c,MP_ROUND)

#define  mp_get_f                      mpfr_get_f

#define  mp_get_str(a,b,c,d,e)         mpfr_get_str(a,b,c,d,e,MP_ROUND)

#define  mp_free_str                   mpfr_free_str

#define  mp_fits_ulong_p(a)            mpfr_fits_ulong_p(a,MP_ROUND)
#define  mp_fits_slong_p(a)            mpfr_fits_slong_p(a,MP_ROUND)
#define  mp_fits_uint_p(a)             mpfr_fits_uint_p(a,MP_ROUND)
#define  mp_fits_sint_p(a)             mpfr_fits_sint_p(a,MP_ROUND)
#define  mp_fits_ushort_p(a)           mpfr_fits_ushort_p(a,MP_ROUND)
#define  mp_fits_sshort_p(a)           mpfr_fits_sshort_p(a,MP_ROUND)
#define  mp_fits_uintmax_p(a)          mpfr_fits_uintmax_p(a,MP_ROUND)
#define  mp_fits_intmax_p(a)           mpfr_fits_intmax_p(a,MP_ROUND)

#define  mp_add(a,b,c)                 mpfr_add(a,b,c,MP_ROUND)
#define  mp_add_ui(a,b,c)              mpfr_add_ui(a,b,c,MP_ROUND)
#define  mp_add_si(a,b,c)              mpfr_add_si(a,b,c,MP_ROUND)
#define  mp_add_d(a,b,c)               mpfr_add_d(a,b,c,MP_ROUND)
#define  mp_add_z(a,b,c)               mpfr_add_z(a,b,c,MP_ROUND)
#define  mp_add_q(a,b,c)               mpfr_add_q(a,b,c,MP_ROUND)

#define  mp_sub(a,b,c)                 mpfr_sub(a,b,c,MP_ROUND)
#define  mp_ui_sub(a,b,c)              mpfr_ui_sub(a,b,c,MP_ROUND)
#define  mp_sub_ui(a,b,c)              mpfr_sub_ui(a,b,c,MP_ROUND)
#define  mp_si_sub(a,b,c)              mpfr_si_sub(a,b,c,MP_ROUND)
#define  mp_sub_si(a,b,c)              mpfr_sub_si(a,b,c,MP_ROUND)
#define  mp_d_sub(a,b,c)               mpfr_d_sub(a,b,c,MP_ROUND)
#define  mp_sub_d(a,b,c)               mpfr_sub_d(a,b,c,MP_ROUND)
#define  mp_z_sub(a,b,c)               mpfr_z_sub(a,b,c,MP_ROUND)
#define  mp_sub_z(a,b,c)               mpfr_sub_z(a,b,c,MP_ROUND)
#define  mp_sub_q(a,b,c)               mpfr_sub_q(a,b,c,MP_ROUND)

#define  mp_mul(a,b,c)                 mpfr_mul(a,b,c,MP_ROUND)
#define  mp_mul_ui(a,b,c)              mpfr_mul_ui(a,b,c,MP_ROUND)
#define  mp_mul_si(a,b,c)              mpfr_mul_si(a,b,c,MP_ROUND)
#define  mp_mul_d(a,b,c)               mpfr_mul_d(a,b,c,MP_ROUND)
#define  mp_mul_z(a,b,c)               mpfr_mul_z(a,b,c,MP_ROUND)
#define  mp_mul_q(a,b,c)               mpfr_mul_q(a,b,c,MP_ROUND)
#define  mp_mul_2exp(a,b,c)            mpfr_mul_2exp(a,b,c,MP_ROUND)
#define  mp_mul_2ui(a,b,c)             mpfr_mul_2ui(a,b,c,MP_ROUND)
#define  mp_mul_2si(a,b,c)             mpfr_mul_2si(a,b,c,MP_ROUND)

#define  mp_sqr(a,b)                   mpfr_sqr(a,b,MP_ROUND)

#define  mp_div(a,b,c)                 mpfr_div(a,b,c,MP_ROUND)
#define  mp_ui_div(a,b,c)              mpfr_ui_div(a,b,c,MP_ROUND)
#define  mp_div_ui(a,b,c)              mpfr_div_ui(a,b,c,MP_ROUND)
#define  mp_si_div(a,b,c)              mpfr_si_div(a,b,c,MP_ROUND)
#define  mp_div_si(a,b,c)              mpfr_div_si(a,b,c,MP_ROUND)
#define  mp_d_div(a,b,c)               mpfr_d_div(a,b,c,MP_ROUND)
#define  mp_div_d(a,b,c)               mpfr_div_d(a,b,c,MP_ROUND)
#define  mp_div_z(a,b,c)               mpfr_div_z(a,b,c,MP_ROUND)
#define  mp_div_q(a,b,c)               mpfr_div_q(a,b,c,MP_ROUND)
#define  mp_div_2exp(a,b,c)            mpfr_div_2exp(a,b,c,MP_ROUND)
#define  mp_div_2ui(a,b,c)             mpfr_div_2ui(a,b,c,MP_ROUND)
#define  mp_div_2si(a,b,c)             mpfr_div_2si(a,b,c,MP_ROUND)

#define  mp_sqrt(a,b)                  mpfr_sqrt(a,b,MP_ROUND)
#define  mp_sqrt_ui(a,b)               mpfr_sqrt_ui(a,b,MP_ROUND)
#define  mp_rec_sqrt(a,b)              mpfr_rec_sqrt(a,b,MP_ROUND)

#define  mp_cbrt(a,b)                  mpfr_cbrt(a,b,MP_ROUND)
#define  mp_root(a,b,c)                mpfr_root(a,b,c,MP_ROUND)

#define  mp_pow(a,b,c)                 mpfr_pow(a,b,c,MP_ROUND)
#define  mp_pow_ui(a,b,c)              mpfr_pow_ui(a,b,c,MP_ROUND)
#define  mp_pow_si(a,b,c)              mpfr_pow_si(a,b,c,MP_ROUND)
#define  mp_pow_z(a,b,c)               mpfr_pow_z(a,b,c,MP_ROUND)
#define  mp_ui_pow_ui(a,b,c)           mpfr_ui_pow_ui(a,b,c,MP_ROUND)
#define  mp_ui_pow(a,b,c)              mpfr_ui_pow(a,b,c,MP_ROUND)

#define  mp_neg(a,b)                   mpfr_neg(a,b,MP_ROUND)
#define  mp_abs(a,b)                   mpfr_abs(a,b,MP_ROUND)

#define  mp_cmp                        mpfr_cmp
#define  mp_cmp_ui                     mpfr_cmp_ui
#define  mp_cmp_si                     mpfr_cmp_si
#define  mp_cmp_d                      mpfr_cmp_d
#define  mp_cmp_ld                     mpfr_cmp_ld
#define  mp_cmp_z                      mpfr_cmp_z
#define  mp_cmp_q                      mpfr_cmp_q
#define  mp_cmp_f                      mpfr_cmp_f

#define  mp_cmp_ui_2exp                mpfr_cmp_ui_2exp
#define  mp_cmp_si_2exp                mpfr_cmp_si_2exp

#define  mp_cmpabs                     mpfr_cmpabs

#define  mp_nan_p                      mpfr_nan_p
#define  mp_inf_p                      mpfr_inf_p
#define  mp_number_p                   mpfr_number_p
#define  mp_zero_p                     mpfr_zero_p
#define  mp_regular_p                  mpfr_regular_p

#define  mp_sgn                        mpfr_sgn

#define  mp_greater_p                  mpfr_greater_p
#define  mp_greaterequal_p             mpfr_greaterequal_p
#define  mp_less_p                     mpfr_less_p
#define  mp_lessequal_p                mpfr_lessequal_p
#define  mp_equal_p                    mpfr_equal_p
#define  mp_lessgreater_p              mpfr_lessgreater_p

#define  mp_unordered_p                mpfr_unordered_p

#define  mp_log(a,b)                   mpfr_log(a,b,MP_ROUND)
#define  mp_log2(a,b)                  mpfr_log2(a,b,MP_ROUND)
#define  mp_log10(a,b)                 mpfr_log10(a,b,MP_ROUND)
#define  mp_log1p                      mpfr_log1p

#define  mp_exp(a,b)                   mpfr_exp(a,b,MP_ROUND)
#define  mp_exp2(a,b)                  mpfr_exp2(a,b,MP_ROUND)
#define  mp_exp10(a,b)                 mpfr_exp10(a,b,MP_ROUND)
#define  mp_expm1(a,b)                 mpfr_expm1(a,b,MP_ROUND)
#define  mp_eint(a,b)                  mpfr_eint(a,b,MP_ROUND)
#define  mp_li2(a,b)                   mpfr_li2(a,b,MP_ROUND)

#define  mp_cos(a,b)                   mpfr_cos(a,b,MP_ROUND)
#define  mp_sin(a,b)                   mpfr_sin(a,b,MP_ROUND)
#define  mp_tan(a,b)                   mpfr_tan(a,b,MP_ROUND)
#define  mp_sin_cos(a,b,c)             mpfr_sin_cos(a,b,c,MP_ROUND)
#define  mp_sec(a,b)                   mpfr_sec(a,b,MP_ROUND)
#define  mp_csc(a,b)                   mpfr_csc(a,b,MP_ROUND)
#define  mp_cot(a,b)                   mpfr_cot(a,b,MP_ROUND)
#define  mp_acos(a,b)                  mpfr_acos(a,b,MP_ROUND)
#define  mp_asin(a,b)                  mpfr_asin(a,b,MP_ROUND)
#define  mp_atan(a,b)                  mpfr_atan(a,b,MP_ROUND)
#define  mp_atan2(a,b,c)               mpfr_atan2(a,b,c,MP_ROUND)
#define  mp_cosh(a,b)                  mpfr_cosh(a,b,MP_ROUND)
#define  mp_sinh(a,b)                  mpfr_sinh(a,b,MP_ROUND)
#define  mp_tanh(a,b)                  mpfr_tanh(a,b,MP_ROUND)
#define  mp_sinh_cosh(a,b,c)           mpfr_sinh_cosh(a,b,c,MP_ROUND)
#define  mp_sech(a,b)                  mpfr_sech(a,b,MP_ROUND)
#define  mp_csch(a,b)                  mpfr_csch(a,b,MP_ROUND)
#define  mp_coth(a,b)                  mpfr_coth(a,b,MP_ROUND)
#define  mp_acosh(a,b)                 mpfr_acosh(a,b,MP_ROUND)
#define  mp_asinh(a,b)                 mpfr_asinh(a,b,MP_ROUND)
#define  mp_atanh(a,b)                 mpfr_atanh(a,b,MP_ROUND)
#define  mp_hypot(a,b,c)               mpfr_hypot(a,b,c,MP_ROUND)

#define  mp_fac_ui                     mpfr_fac_ui

#define  mp_gamma(a,b)                 mpfr_gamma(a,b,MP_ROUND)
#define  mp_lngamma(a,b)               mpfr_lngamma(a,b,MP_ROUND)
#define  mp_lgamma(a,b,c)              mpfr_lgamma(a,b,c,MP_ROUND)
#define  mp_digamma(a,b)               mpfr_digamma(a,b,MP_ROUND)

#define  mp_zeta(a,b)                  mpfr_zeta(a,b,MP_ROUND)
#define  mp_zeta_ui(a,b)               mpfr_zeta_ui(a,b,MP_ROUND)

#define  mp_erf(a,b)                   mpfr_erf(a,b,MP_ROUND)
#define  mp_erfc(a,b)                  mpfr_erfc(a,b,MP_ROUND)

#define  mp_j0(a,b)                    mpfr_j0(a,b,MP_ROUND)
#define  mp_j1(a,b)                    mpfr_j1(a,b,MP_ROUND)
#define  mp_jn(a,b,c)                  mpfr_jn(a,b,c,MP_ROUND)

#define  mp_y0(a,b)                    mpfr_y0(a,b,MP_ROUND)
#define  mp_y1(a,b)                    mpfr_y1(a,b,MP_ROUND)
#define  mp_yn(a,b,c)                  mpfr_yn(a,b,c,MP_ROUND)

#define  mp_fma(a,b,c,d)               mpfr_fma(a,b,c,d,MP_ROUND)
#define  mp_fms(a,b,c,d)               mpfr_fms(a,b,c,d,MP_ROUND)

#define  mp_agm(a,b,c)                 mpfr_agm(a,b,c,MP_ROUND)

#define  mp_ai(a,b)                    mpfr_ai(a,b,MP_ROUND)

#define  mp_const_log2(a)              mpfr_const_log2(a,MP_ROUND)
#define  mp_const_pi(a)                mpfr_const_pi(a,MP_ROUND)
#define  mp_const_euler(a)             mpfr_const_euler(a,MP_ROUND)
#define  mp_const_catalan(a)           mpfr_const_catalan(a,MP_ROUND)

#define  mp_free_cache                 mpfr_free_cache

#define  mp_sum(a,b,c)                 mpfr_sum(a,b,c,MP_ROUND)

#define  mp_out_str(a,b,c)             mpfr_out_str(a,b,c,MP_ROUND)
#define  mp_inp_str(a,b,c)             mpfr_inp_str(a,b,c,MP_ROUND)

#define  mp_fprintf                    mpfr_fprintf
#define  mp_vfprintf                   mpfr_vfprintf
#define  mp_printf                     mpfr_printf
#define  mp_vprintf                    mpfr_vprintf
#define  mp_sprintf                    mpfr_sprintf
#define  mp_vsprintf                   mpfr_vsprintf
#define  mp_snprintf                   mpfr_snprintf
#define  mp_vsnprintf                  mpfr_vsnprintf
#define  mp_asprintf                   mpfr_asprintf
#define  mp_vasprintf                  mpfr_vasprintf

#define  mp_ceil                       mpfr_ceil
#define  mp_floor                      mpfr_floor
#define  mp_round                      mpfr_round
#define  mp_trunc                      mpfr_trunc

#define  mp_rint(a,b)                  mpfr_rint(a,b,MP_ROUND)
#define  mp_rint_ceil(a,b)             mpfr_rint_ceil(a,b,MP_ROUND)
#define  mp_rint_floor(a,b)            mpfr_rint_floor(a,b,MP_ROUND)
#define  mp_rint_round(a,b)            mpfr_rint_round(a,b,MP_ROUND)
#define  mp_rint_trunc(a,b)            mpfr_rint_trunc(a,b,MP_ROUND)

#define  mp_frac(a,b)                  mpfr_frac(a,b,MP_ROUND)
#define  mp_modf(a,b,c)                mpfr_modf(a,b,c,MP_ROUND)
#define  mp_fmod(a,b,c)                mpfr_fmod(a,b,c,MP_ROUND)
#define  mp_remainder(a,b,c)           mpfr_remainder(a,b,c,MP_ROUND)
#define  mp_remquo(a,b,c,d)            mpfr_remquo(a,b,c,d,MP_ROUND)

#define  mp_integer_p                  mpfr_integer_p

#define  mp_set_default_rounding_mode  mpfr_set_default_rounding_mode
#define  mp_get_default_rounding_mode  mpfr_get_default_rounding_mode
#define  mp_prec_round                 mpfr_prec_round
#define  mp_can_round                  mpfr_can_round
#define  mp_min_prec                   mpfr_min_prec
#define  mp_print_rnd_mode             mpfr_print_rnd_mode

#define  mp_nexttoward                 mpfr_nexttoward
#define  mp_nextabove                  mpfr_nextabove
#define  mp_nextbelow                  mpfr_nextbelow

#define  mp_min(a,b,c)                 mpfr_min(a,b,c,MP_ROUND)
#define  mp_max(a,b,c)                 mpfr_max(a,b,c,MP_ROUND)
#define  mp_dim(a,b,c)                 mpfr_dim(a,b,c,MP_ROUND)

#define  mp_urandomb(a)                mpfr_urandomb(a,MP_ROUND)
#define  mp_urandom(a,b)               mpfr_urandom(a,b,MP_ROUND)
#define  mp_grandom(a,b,c)             mpfr_grandom(a,b,c,MP_ROUND)

#define  mp_get_exp                    mpfr_get_exp
#define  mp_set_exp                    mpfr_set_exp

#define  mp_signbit                    mpfr_signbit
#define  mp_setsign(a,b)               mpfr_setsign(a,b,MP_ROUND)
#define  mp_copysign(a,b,c)            mpfr_copysign(a,b,c,MP_ROUND)

#define  mp_get_version                mpfr_get_version
#define  mp_get_patches                mpfr_get_patches
#define  mp_buildopt_tls_p             mpfr_buildopt_tls_p
#define  mp_buildopt_decimal_p         mpfr_buildopt_decimal_p
#define  mp_buildopt_gmpinternals_p    mpfr_buildopt_gmpinternals_p
#define  mp_buildopt_tune_case         mpfr_buildopt_tune_case

#define  mp_get_emin                   mpfr_get_emin
#define  mp_get_emax                   mpfr_get_emax
#define  mp_set_emin                   mpfr_set_emin
#define  mp_set_emax                   mpfr_set_emax
#define  mp_get_emin_min               mpfr_get_emin_min
#define  mp_get_emin_max               mpfr_get_emin_max
#define  mp_get_emax_min               mpfr_get_emax_min
#define  mp_get_emax_max               mpfr_get_emax_max

#define  mp_check_range                mpfr_check_range
#define  mp_subnormalize(a,b)          mpfr_subnormalize(a,b,MP_ROUND)

#define  mp_clear_underflow            mpfr_clear_underflow
#define  mp_clear_overflow             mpfr_clear_overflow
#define  mp_clear_divby0               mpfr_clear_divby0
#define  mp_clear_nanflag              mpfr_clear_nanflag
#define  mp_clear_inexflag             mpfr_clear_inexflag
#define  mp_clear_erangeflag           mpfr_clear_erangeflag
#define  mp_set_underflow              mpfr_set_underflow
#define  mp_set_overflow               mpfr_set_overflow
#define  mp_set_divby0                 mpfr_set_divby0
#define  mp_set_nanflag                mpfr_set_nanflag
#define  mp_set_inexflag               mpfr_set_inexflag
#define  mp_set_erangeflag             mpfr_set_erangeflag
#define  mp_clear_flags                mpfr_clear_flags
#define  mp_underflow_p                mpfr_underflow_p
#define  mp_overflow_p                 mpfr_overflow_p
#define  mp_divby0_p                   mpfr_divby0_p
#define  mp_nanflag_p                  mpfr_nanflag_p
#define  mp_inexflag_p                 mpfr_inexflag_p
#define  mp_erangeflag_p               mpfr_erangeflag_p

#else

// If not actually using MPFR, then emulate it by direct floating-point
// arithmetic statements.

typedef  long double  mp_real;

#define  mp_init2(a,b)           assert((b) <= 80), (a) = 0
#define  mp_clear(a)             (a) = 0
#define  mp_get_prec(a)          64
#define  mp_set(a,b)             (a) = (b);
#define  mp_set_d(a,b)           (a) = (mp_real)(b)
#define  mp_set_str(a,b,c)       assert((c) == 10), (a) = strtold((b), NULL)
#define  mp_get_d(a)             (double)(a)
#define  mp_add(a,b,c)           (a) = (b) + (c)
#define  mp_add_d(a,b,c)         (a) = (b) + (mp_real)(c)
#define  mp_sub(a,b,c)           (a) = (b) - (c)
#define  mp_sub_d(a,b,c)         (a) = (b) - (mp_real)(c)
#define  mp_mul(a,b,c)           (a) = (b) * (c)
#define  mp_mul_d(a,b,c)         (a) = (b) * (mp_real)(c)
#define  mp_div(a,b,c)           (a) = (b) / (c)
#define  mp_div_d(a,b,c)         (a) = (b) / (mp_real)(c)
#define  mp_sqr(a,b)             (a) = (b) * (b)
#define  mp_sqrt(a,b)            (a) = sqrtl(b);
#define  mp_log2(a,b)            (a) = log2l(b);
#define  mp_abs(a,b)             (a) = fabsl(b);
#define  mp_sgn(a)               (((a) > 0) - ((a) < 0))
#define  mp_nan_p(a)             isnan(a)
#define  mp_inf_p(a)             isinf(a)
#define  mp_zero_p(a)            ((a) == 0)
#define  mp_equal_p(a,b)         ((a) == (b))
#define  mp_greater_p(a,b)       ((a) >  (b))
#define  mp_greaterequal_p(a,b)  ((a) >= (b))
#define  mp_less_p(a,b)          ((a) <  (b))
#define  mp_lessequal_p(a,b)     ((a) <= (b))

//#define  mp_fprintf              ()
//#define  mp_printf               ()
//#define  mp_snprintf             ()

#endif


//-----------------------------------------------------------------------------
// EXTENSIONS

extern_public_function
  void mp_lerp(mp_real *x,
               const mp_real t, const mp_real x0, const mp_real x1);

extern_public_function
  void mp_lerp_d(mp_real *x,
                 const real t, const mp_real x0, const mp_real x1);

extern_public_function
  int mp_get_min_prec_from_string(const char *str);
