/*################################################################################
  ##
  ##   Copyright (C) 2011-2023 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * pdf of the Poisson distribution
 */

//
// scalar input

namespace internal
{

template<typename T>
statslib_constexpr
T
dpois_log_compute(const llint_t x, const T rate_par)
noexcept
{
    return( x * stmath::log(rate_par) - rate_par - stmath::lgamma(T(x+1)) );
}

template<typename T>
statslib_constexpr
T
dpois_vals_check(const llint_t x, const T rate_par, const bool log_form)
noexcept
{
    return( !pois_sanity_check(rate_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            x < llint_t(0) ? \
                log_zero_if<T>(log_form) :
            //
            rate_par == T(0) ? \
                x == llint_t(0) ? \
                    log_one_if<T>(log_form) :
                    log_zero_if<T>(log_form) :
            //
            GCINT::is_posinf(rate_par) ? \
                log_zero_if<T>(log_form) :
            //
            exp_if(dpois_log_compute(x,rate_par), !log_form) );
}

}

template<typename T>
statslib_constexpr
return_t<T>
dpois(const llint_t x, const T rate_par, const bool log_form)
noexcept
{
    return internal::dpois_vals_check(x,static_cast<return_t<T>>(rate_par),log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
dpois_vec(const eT* __stats_pointer_settings__ vals_in, const T1 rate_par, const bool log_form, 
                rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dpois,vals_in,vals_out,num_elem,rate_par,log_form);
}
#endif

}

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
dpois(const std::vector<eT>& x, const T1 rate_par, const bool log_form)
{
    STDVEC_DIST_FN(dpois_vec,rate_par,log_form);
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
dpois(const ArmaMat<eT>& X, const T1 rate_par, const bool log_form)
{
    ARMA_DIST_FN(dpois_vec,rate_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
dpois(const ArmaGen<mT,tT>& X, const T1 rate_par, const bool log_form)
{
    return dpois(X.eval(),rate_par,log_form);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dpois(const BlazeMat<eT,To>& X, const T1 rate_par, const bool log_form)
{
    BLAZE_DIST_FN(dpois,rate_par,log_form);
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dpois(const EigenMat<eT,iTr,iTc>& X, const T1 rate_par, const bool log_form)
{
    EIGEN_DIST_FN(dpois_vec,rate_par,log_form);
}
#endif
