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
 * pdf of the Rademacher distribution
 */

#ifndef _statslib_dradem_HPP
#define _statslib_dradem_HPP

//
// scalar input

/**
 * @brief Density function of the Rademacher distribution
 *
 * @param x an integral-valued input, equal to -1 or 1.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return the density function evaluated at \c x.
 * 
 * Example:
 * \code{.cpp} stats::dradem(1,0.6,false); \endcode
 */

template<typename T>
statslib_constexpr
return_t<T>
dradem(const llint_t x, const T prob_par, const bool log_form = false) noexcept;

//
// vector/matrix input

/**
 * @brief Density function of the Rademacher distribution
 *
 * @param x a standard vector.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a vector of density function values corresponding to the elements of \c x.
 * 
 * Example:
 * \code{.cpp}
 * std::vector<int> x = {-1, 1, 0};
 * stats::dradem(x,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_STDVEC_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>>
statslib_inline
std::vector<rT>
dradem(const std::vector<eT>& x, const T1 prob_par, const bool log_form = false);
#endif

/**
 * @brief Density function of the Rademacher distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {-1, 0,  1},
 *                 {-1, 1, -1} };
 * stats::dradem(X,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>>
statslib_inline
ArmaMat<rT>
dradem(const ArmaMat<eT>& X, const T1 prob_par, const bool log_form = false);

/**
 * @brief Density function of the Rademacher distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * arma::mat X = { {-1, 0,  1},
 *                 {-1, 1, -1} };
 * stats::dradem(X,0.5,false);
 * \endcode
 */

template<typename mT, typename tT, typename T1>
statslib_inline
mT 
dradem(const ArmaGen<mT,tT>& X, const T1 prob_par, const bool log_form = false);
#endif

/**
 * @brief Density function of the Rademacher distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dradem(X,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>, bool To = blaze::columnMajor>
statslib_inline
BlazeMat<rT,To>
dradem(const BlazeMat<eT,To>& X, const T1 prob_par, const bool log_form = false);
#endif

/**
 * @brief Density function of the Rademacher distribution
 *
 * @param X a matrix of input values.
 * @param prob_par the probability parameter, a real-valued input.
 * @param log_form return the log-density or the true form.
 *
 * @return a matrix of density function values corresponding to the elements of \c X.
 * 
 * Example:
 * \code{.cpp}
 * stats::dradem(X,0.5,false);
 * \endcode
 */

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, typename T1, typename rT = common_return_t<eT,T1>, int iTr = Eigen::Dynamic, int iTc = Eigen::Dynamic>
statslib_inline
EigenMat<rT,iTr,iTc>
dradem(const EigenMat<eT,iTr,iTc>& X, const T1 prob_par, const bool log_form = false);
#endif

//
// include implementation files

#include "dradem.ipp"

#endif
