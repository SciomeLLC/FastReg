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
 * cdf of the chi-squared distribution
 */

//
// scalar input

namespace internal
{

template<typename T>
statslib_constexpr
T
pchisq_compute(const T x, const T dof_par)
noexcept
{
    return gcem::incomplete_gamma(dof_par/T(2),x/T(2));
}

template<typename T>
statslib_constexpr
T
pchisq_vals_check(const T x, const T dof_par, const bool log_form)
noexcept
{
    return( !chisq_sanity_check(x,dof_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            STLIM<T>::epsilon() > x ? \
                log_zero_if<T>(log_form) : 
            // now x > 0 cases
            dof_par == T(0) ? \
                log_one_if<T>(log_form) :
            //
            GCINT::is_posinf(x) ? \
                log_one_if<T>(log_form) :
            GCINT::is_posinf(dof_par) ? \
                log_zero_if<T>(log_form) :
            //
            log_if(pchisq_compute(x,dof_par), log_form) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
pchisq_type_check(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return pchisq_vals_check(static_cast<TC>(x),static_cast<TC>(dof_par),log_form);
}

}

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
pchisq(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return internal::pchisq_type_check(x,dof_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
pchisq_vec(const eT* __stats_pointer_settings__ vals_in, const T1 dof_par, const bool log_form, 
                 rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(pchisq,vals_in,vals_out,num_elem,dof_par,log_form);
}
#endif

}

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
pchisq(const std::vector<eT>& x, const T1 dof_par, const bool log_form)
{
    STDVEC_DIST_FN(pchisq_vec,dof_par,log_form);
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
pchisq(const ArmaMat<eT>& X, const T1 dof_par, const bool log_form)
{
    ARMA_DIST_FN(pchisq_vec,dof_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
pchisq(const ArmaGen<mT,tT>& X, const T1 dof_par, const bool log_form)
{
    return pchisq(X.eval(),dof_par,log_form);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
pchisq(const BlazeMat<eT,To>& X, const T1 dof_par, const bool log_form)
{
    BLAZE_DIST_FN(pchisq_vec,dof_par,log_form);
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
pchisq(const EigenMat<eT,iTr,iTc>& X, const T1 dof_par, const bool log_form)
{
    EIGEN_DIST_FN(pchisq_vec,dof_par,log_form);
}
#endif
