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
 * pdf of the Cauchy distribution
 */

//
// scalar input

namespace internal
{

template<typename T>
statslib_constexpr
T
dcauchy_compute(const T z, const T sigma_par)
noexcept
{
    return( T(1) / ( sigma_par*GCEM_PI*(T(1) + z*z) ) );
}

template<typename T>
statslib_constexpr
T
dcauchy_limit_vals(const T x, const T mu_par, const T sigma_par)
noexcept
{
    return( // sigma == Inf
            GCINT::is_posinf(sigma_par) ? \
                GCINT::any_inf(x,mu_par) ? \
                    STLIM<T>::quiet_NaN() :
                    T(0) :
            // sigma finite; x == mu == Inf or -Inf 
            GCINT::all_posinf(x,mu_par) || GCINT::all_neginf(x,mu_par) ? \
                STLIM<T>::quiet_NaN() :
            // sigma finite; x and mu have opposite Inf signs
            GCINT::is_posinf(x) && GCINT::is_neginf(mu_par) ? \
                T(0) :
            GCINT::is_neginf(x) && GCINT::is_posinf(mu_par) ? \
                T(0) :
            //
                T(0) );
}

template<typename T>
statslib_constexpr
T
dcauchy_vals_check(const T x, const T mu_par, const T sigma_par, const bool log_form)
noexcept
{
    return( !cauchy_sanity_check(x,mu_par,sigma_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::any_inf(x,mu_par,sigma_par) ? \
                log_if(dcauchy_limit_vals(x,mu_par,sigma_par),log_form) :
            //
            log_if(dcauchy_compute((x-mu_par)/sigma_par,sigma_par), log_form) );
}

template<typename T1, typename T2, typename T3, typename TC = common_return_t<T1,T2,T3>>
statslib_constexpr
TC
dcauchy_type_check(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return dcauchy_vals_check(static_cast<TC>(x),static_cast<TC>(mu_par),
                              static_cast<TC>(sigma_par),log_form);
}

}

template<typename T1, typename T2, typename T3>
statslib_constexpr
common_return_t<T1,T2,T3>
dcauchy(const T1 x, const T2 mu_par, const T3 sigma_par, const bool log_form)
noexcept
{
    return internal::dcauchy_type_check(x,mu_par,sigma_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
void
dcauchy_vec(const eT* __stats_pointer_settings__ vals_in, const T1 mu_par, const T2 sigma_par, const bool log_form, 
                  rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dcauchy,vals_in,vals_out,num_elem,mu_par,sigma_par,log_form);
}
#endif

}

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
std::vector<rT>
dcauchy(const std::vector<eT>& x, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    STDVEC_DIST_FN(dcauchy_vec,mu_par,sigma_par,log_form);
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT>
statslib_inline
ArmaMat<rT>
dcauchy(const ArmaMat<eT>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    ARMA_DIST_FN(dcauchy_vec,mu_par,sigma_par,log_form);
}

template<typename mT, typename tT, typename Tb>
statslib_inline
mT
dcauchy(const ArmaGen<mT,tT>& X, const Tb mu_par, const Tb sigma_par, const bool log_form)
{
    return dcauchy(X.eval(),mu_par,sigma_par,log_form);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dcauchy(const BlazeMat<eT,To>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    BLAZE_DIST_FN(dcauchy,mu_par,sigma_par,log_form);
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename T2, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dcauchy(const EigenMat<eT,iTr,iTc>& X, const T1 mu_par, const T2 sigma_par, const bool log_form)
{
    EIGEN_DIST_FN(dcauchy_vec,mu_par,sigma_par,log_form);
}
#endif
