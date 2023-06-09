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
 * cdf of the univariate logistic distribution
 */

//
// scalar input

namespace internal
{

template<typename T>
statslib_constexpr
T
plogis_compute(const T z)
{
    return( T(0.5)*( T(1) + gcem::tanh(z) ) );
}

template<typename T>
statslib_constexpr
T
plogis_limit_vals(const T x, const T mu_par, const T sigma_par)
noexcept
{
    return( // x == mu == Inf or -Inf 
            GCINT::all_posinf(x,mu_par) || GCINT::all_neginf(x,mu_par) ? \
                STLIM<T>::quiet_NaN() :
            // sigma == Inf
            GCINT::is_posinf(sigma_par) ? \
                // x or mu == +/-Inf
                GCINT::any_inf(x,mu_par) ? \
                    STLIM<T>::quiet_NaN() :
                // x and mu finite
                    T(0.5) :
            // sigma finite; x and mu can have opposite Inf signs
            GCINT::is_posinf(x) || GCINT::is_neginf(mu_par) ? \
                T(1) :
            GCINT::is_neginf(x) || GCINT::is_posinf(mu_par) ? \
                T(0) :
            // sigma == 0
            sigma_par == T(0) ? \
                T(0.5) + T(0.5)*gcem::sgn(x-mu_par) : // arbitrary convention?
            //
                T(0) );
}

template<typename T>
statslib_constexpr
T
plogis_vals_check(const T x, const T mu_par, const T sigma_par, const bool log_form)
{
    return( !logis_sanity_check(x,mu_par,sigma_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::any_inf(x,mu_par,sigma_par) || sigma_par == T(0) ? \
                log_if(plogis_limit_vals(x,mu_par,sigma_par),log_form) :
            //
            log_if(plogis_compute((x-mu_par)/(2*sigma_par)), log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
plogis_type_check(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return plogis_vals_check(static_cast<TC>(x),static_cast<TC>(mu_par),
                             static_cast<TC>(sigma_par),log_form);
}

}

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
plogis(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return internal::plogis_type_check(x,mu_par,sigma_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
plogis_vec(const eT* __stats_pointer_settings__ vals_in, const T1 mu_par, const T2 sigma_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(plogis,vals_in,vals_out,num_elem,mu_par,sigma_par,log_form);
}
#endif

}

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
plogis(const std::vector<eT>& x, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    STDVEC_DIST_FN(plogis_vec,mu_par,sigma_par,log_form);
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
plogis(const ArmaMat<eT>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    ARMA_DIST_FN(plogis_vec,mu_par,sigma_par,log_form);
}

template<typename mT, typename tT, typename T1, typename T2>
statslib_inline
mT
plogis(const ArmaGen<mT,tT>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    return plogis(X.eval(),mu_par,sigma_par,log_form);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
plogis(const BlazeMat<eT,To>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    BLAZE_DIST_FN(plogis_vec,mu_par,sigma_par,log_form);
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
plogis(const EigenMat<eT,iTr,iTc>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    EIGEN_DIST_FN(plogis_vec,mu_par,sigma_par,log_form);
}
#endif
