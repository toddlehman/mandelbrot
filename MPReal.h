/*-----------------------------------------------------------------------------
  Copyright (c) 2013 by Todd S. Lehman.  All rights reserved.
-----------------------------------------------------------------------------*/

#import "Common.h"


//------------------------------------------------------------------------------
// Include MPFR multiprecision floating-point library and modify its interface
// to be a bit cleaner.

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

#define  mp_set                        mpfr_set                      
#define  mp_set_ui                     mpfr_set_ui                   
#define  mp_set_si                     mpfr_set_si                   
#define  mp_set_uj                     mpfr_set_uj                   
#define  mp_set_sj                     mpfr_set_sj                   
#define  mp_set_flt                    mpfr_set_flt                  
#define  mp_set_d                      mpfr_set_d                    
#define  mp_set_ld                     mpfr_set_ld                   
#define  mp_set_decimal64              mpfr_set_decimal64            
#define  mp_set_z                      mpfr_set_z                    
#define  mp_set_q                      mpfr_set_q                    
#define  mp_set_f                      mpfr_set_f                    

#define  mp_set_ui_2exp                mpfr_set_ui_2exp              
#define  mp_set_si_2exp                mpfr_set_si_2exp              
#define  mp_set_uj_2exp                mpfr_set_uj_2exp              
#define  mp_set_sj_2exp                mpfr_set_sj_2exp              
#define  mp_set_z_2exp                 mpfr_set_z_2exp               

#define  mp_set_str                    mpfr_set_str                  

#define  mp_strtofr                    mpfr_strtofr                  

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

#define  mp_get_flt                    mpfr_get_flt                  
#define  mp_get_d                      mpfr_get_d                    
#define  mp_get_ld                     mpfr_get_ld                   
#define  mp_get_decimal64              mpfr_get_decimal64            

#define  mp_get_si                     mpfr_get_si                   
#define  mp_get_ui                     mpfr_get_ui                   
#define  mp_get_sj                     mpfr_get_sj                   
#define  mp_get_uj                     mpfr_get_uj                   

#define  mp_get_d_2exp                 mpfr_get_d_2exp               
#define  mp_get_ld_2exp                mpfr_get_ld_2exp              

#define  mp_frexp                      mpfr_frexp                    

#define  mp_get_z_2exp                 mpfr_get_z_2exp               

#define  mp_get_z                      mpfr_get_z                    

#define  mp_get_f                      mpfr_get_f                    

#define  mp_get_str                    mpfr_get_str                  

#define  mp_free_str                   mpfr_free_str                 

#define  mp_fits_ulong_p               mpfr_fits_ulong_p             
#define  mp_fits_slong_p               mpfr_fits_slong_p             
#define  mp_fits_uint_p                mpfr_fits_uint_p              
#define  mp_fits_sint_p                mpfr_fits_sint_p              
#define  mp_fits_ushort_p              mpfr_fits_ushort_p            
#define  mp_fits_sshort_p              mpfr_fits_sshort_p            
#define  mp_fits_uintmax_p             mpfr_fits_uintmax_p           
#define  mp_fits_intmax_p              mpfr_fits_intmax_p            

#define  mp_add                        mpfr_add                      
#define  mp_add_ui                     mpfr_add_ui                   
#define  mp_add_si                     mpfr_add_si                   
#define  mp_add_d                      mpfr_add_d                    
#define  mp_add_z                      mpfr_add_z                    
#define  mp_add_q                      mpfr_add_q                    

#define  mp_sub                        mpfr_sub                      
#define  mp_ui_sub                     mpfr_ui_sub                   
#define  mp_sub_ui                     mpfr_sub_ui                   
#define  mp_si_sub                     mpfr_si_sub                   
#define  mp_sub_si                     mpfr_sub_si                   
#define  mp_d_sub                      mpfr_d_sub                    
#define  mp_sub_d                      mpfr_sub_d                    
#define  mp_z_sub                      mpfr_z_sub                    
#define  mp_sub_z                      mpfr_sub_z                    
#define  mp_sub_q                      mpfr_sub_q                    

#define  mp_mul                        mpfr_mul                      
#define  mp_mul_ui                     mpfr_mul_ui                   
#define  mp_mul_si                     mpfr_mul_si                   
#define  mp_mul_d                      mpfr_mul_d                    
#define  mp_mul_z                      mpfr_mul_z                    
#define  mp_mul_q                      mpfr_mul_q                    

#define  mp_sqr                        mpfr_sqr                      

#define  mp_div                        mpfr_div                      
#define  mp_ui_div                     mpfr_ui_div                   
#define  mp_div_ui                     mpfr_div_ui                   
#define  mp_si_div                     mpfr_si_div                   
#define  mp_div_si                     mpfr_div_si                   
#define  mp_d_div                      mpfr_d_div                    
#define  mp_div_d                      mpfr_div_d                    
#define  mp_div_z                      mpfr_div_z                    
#define  mp_div_q                      mpfr_div_q                    

#define  mp_sqrt                       mpfr_sqrt                     
#define  mp_sqrt_ui                    mpfr_sqrt_ui                  

#define  mp_rec_sqrt                   mpfr_rec_sqrt                 

#define  mp_cbrt                       mpfr_cbrt                     
#define  mp_root                       mpfr_root                     

#define  mp_pow                        mpfr_pow                      
#define  mp_pow_ui                     mpfr_pow_ui                   
#define  mp_pow_si                     mpfr_pow_si                   
#define  mp_pow_z                      mpfr_pow_z                    
#define  mp_ui_pow_ui                  mpfr_ui_pow_ui                
#define  mp_ui_pow                     mpfr_ui_pow                   

#define  mp_neg                        mpfr_neg                      
#define  mp_abs                        mpfr_abs                      

#define  mp_dim                        mpfr_dim                      

#define  mp_mul_2ui                    mpfr_mul_2ui                  
#define  mp_mul_2si                    mpfr_mul_2si                  

#define  mp_div_2ui                    mpfr_div_2ui                  
#define  mp_div_2si                    mpfr_div_2si                  

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

#define  mp_log                        mpfr_log                      
#define  mp_log2                       mpfr_log2                     
#define  mp_log10                      mpfr_log10                    

#define  mp_exp                        mpfr_exp                      
#define  mp_exp2                       mpfr_exp2                     
#define  mp_exp10                      mpfr_exp10                    

#define  mp_cos                        mpfr_cos                      
#define  mp_sin                        mpfr_sin                      
#define  mp_tan                        mpfr_tan                      
#define  mp_sin_cos                    mpfr_sin_cos                  
#define  mp_sec                        mpfr_sec                      
#define  mp_csc                        mpfr_csc                      
#define  mp_cot                        mpfr_cot                      
#define  mp_acos                       mpfr_acos                     
#define  mp_asin                       mpfr_asin                     
#define  mp_atan                       mpfr_atan                     
#define  mp_atan2                      mpfr_atan2                    
#define  mp_cosh                       mpfr_cosh                     
#define  mp_sinh                       mpfr_sinh                     
#define  mp_tanh                       mpfr_tanh                     
#define  mp_sinh_cosh                  mpfr_sinh_cosh                
#define  mp_sech                       mpfr_sech                     
#define  mp_csch                       mpfr_csch                     
#define  mp_coth                       mpfr_coth                     
#define  mp_acosh                      mpfr_acosh                    
#define  mp_asinh                      mpfr_asinh                    
#define  mp_atanh                      mpfr_atanh                    

#define  mp_fac_ui                     mpfr_fac_ui                   

#define  mp_log1p                      mpfr_log1p                    

#define  mp_expm1                      mpfr_expm1                    

#define  mp_eint                       mpfr_eint                     

#define  mp_li2                        mpfr_li2                      

#define  mp_gamma                      mpfr_gamma                    
#define  mp_lngamma                    mpfr_lngamma                  
#define  mp_lgamma                     mpfr_lgamma                   
#define  mp_digamma                    mpfr_digamma                  

#define  mp_zeta                       mpfr_zeta                     
#define  mp_zeta_ui                    mpfr_zeta_ui                  

#define  mp_erf                        mpfr_erf                      
#define  mp_erfc                       mpfr_erfc                     

#define  mp_j0                         mpfr_j0                       
#define  mp_j1                         mpfr_j1                       
#define  mp_jn                         mpfr_jn                       

#define  mp_y0                         mpfr_y0                       
#define  mp_y1                         mpfr_y1                       
#define  mp_yn                         mpfr_yn                       

#define  mp_fma                        mpfr_fma                      
#define  mp_fms                        mpfr_fms                      

#define  mp_agm                        mpfr_agm                      

#define  mp_hypot                      mpfr_hypot                    

#define  mp_ai                         mpfr_ai                       

#define  mp_const_log2                 mpfr_const_log2               
#define  mp_const_pi                   mpfr_const_pi                 
#define  mp_const_euler                mpfr_const_euler              
#define  mp_const_catalan              mpfr_const_catalan            

#define  mp_free_cache                 mpfr_free_cache               

#define  mp_sum                        mpfr_sum                      

#define  mp_out_str                    mpfr_out_str                  
#define  mp_inp_str                    mpfr_inp_str                  

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

#define  mp_rint                       mpfr_rint                     
#define  mp_ceil                       mpfr_ceil                     
#define  mp_floor                      mpfr_floor                    
#define  mp_round                      mpfr_round                    
#define  mp_trunc                      mpfr_trunc                    

#define  mp_rint_ceil                  mpfr_rint_ceil                
#define  mp_rint_floor                 mpfr_rint_floor               
#define  mp_rint_round                 mpfr_rint_round               
#define  mp_rint_trunc                 mpfr_rint_trunc               

#define  mp_frac                       mpfr_frac                     
#define  mp_modf                       mpfr_modf                     
#define  mp_fmod                       mpfr_fmod                     
#define  mp_remainder                  mpfr_remainder                
#define  mp_remquo                     mpfr_remquo                   

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

#define  mp_min                        mpfr_min                      
#define  mp_max                        mpfr_max                      

#define  mp_urandomb                   mpfr_urandomb                 
#define  mp_urandom                    mpfr_urandom                  
#define  mp_grandom                    mpfr_grandom                  

#define  mp_get_exp                    mpfr_get_exp                  
#define  mp_set_exp                    mpfr_set_exp                  

#define  mp_signbit                    mpfr_signbit                  
#define  mp_setsign                    mpfr_setsign                  
#define  mp_copysign                   mpfr_copysign                 

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
#define  mp_subnormalize               mpfr_subnormalize             

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
