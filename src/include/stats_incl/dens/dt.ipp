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
 * pdf of the t-distribution
 */

//
// scalar input

namespace internal
{

template<typename T>
statslib_constexpr
T
dt_log_mult_term(const T z, const T dof_par)
noexcept
{
    return( - (dof_par/T(2) + T(0.5)) * stmath::log(T(1) + (z/dof_par)*z) );
}

template<typename T>
statslib_constexpr
T
dt_log_cons_term(const T dof_par)
noexcept
{
    return( stmath::lgamma(dof_par/T(2) + T(0.5)) \
                - T(0.5)*( stmath::log(dof_par) + T(GCEM_LOG_PI) ) \
                - stmath::lgamma(dof_par/T(2)) );
}

template<typename T>
statslib_constexpr
T
dt_log_compute(const T z, const T dof_par)
noexcept
{
    return( dt_log_cons_term(dof_par) + dt_log_mult_term(z,dof_par) );
}

template<typename T>
statslib_constexpr
T
dt_vals_check(const T x, const T dof_par, const bool log_form)
noexcept
{
    return( !t_sanity_check(x,dof_par) ? \
                STLIM<T>::quiet_NaN() :
            //
            GCINT::is_inf(x) ? \
                log_zero_if<T>(log_form) :
            //
            GCINT::is_posinf(dof_par) ? \
                dnorm(x,T(0),T(1),log_form) :
            //
            exp_if(dt_log_compute(x,dof_par), !log_form) );
}

template<typename T1, typename T2, typename TC = common_return_t<T1,T2>>
statslib_constexpr
TC
dt_type_check(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return dt_vals_check(static_cast<TC>(x),static_cast<TC>(dof_par),log_form);
}

}

template<typename T1, typename T2>
statslib_constexpr
common_return_t<T1,T2>
dt(const T1 x, const T2 dof_par, const bool log_form)
noexcept
{
    return internal::dt_type_check(x,dof_par,log_form);
}

//
// vector/matrix input

namespace internal
{

#ifdef STATS_ENABLE_INTERNAL_VEC_FEATURES
template<typename eT, typename T1, typename rT>
statslib_inline
void
dt_vec(const eT* __stats_pointer_settings__ vals_in, const T1 dof_par, const bool log_form, 
             rT* __stats_pointer_settings__ vals_out, const ullint_t num_elem)
{
    EVAL_DIST_FN_VEC(dt,vals_in,vals_out,num_elem,dof_par,log_form);
}
#endif

}

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
std::vector<rT>
dt(const std::vector<eT>& x, const T1 dof_par, const bool log_form)
{
    STDVEC_DIST_FN(dt_vec,dof_par,log_form);
}
#endif

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT>
statslib_inline
ArmaMat<rT>
dt(const ArmaMat<eT>& X, const T1 dof_par, const bool log_form)
{
    ARMA_DIST_FN(dt_vec,dof_par,log_form);
}

template<typename mT, typename tT, typename T1>
statslib_inline
mT
dt(const ArmaGen<mT,tT>& X, const T1 dof_par, const bool log_form)
{
    return dt(X.eval(),dof_par,log_form);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT, bool To>
statslib_inline
BlazeMat<rT,To>
dt(const BlazeMat<eT,To>& X, const T1 dof_par, const bool log_form)
{
    BLAZE_DIST_FN(dt,dof_par,log_form);
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT, int iTr, int iTc>
statslib_inline
EigenMat<rT,iTr,iTc>
dt(const EigenMat<eT,iTr,iTc>& X, const T1 dof_par, const bool log_form)
{
    EIGEN_DIST_FN(dt_vec,dof_par,log_form);
}
#endif
