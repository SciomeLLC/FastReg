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
 * for internal use only; used to switch between the different matrix libraries
 */

//
// identity matrix

#ifdef STATS_ENABLE_ARMA_WRAPPERS
template<typename eT>
statslib_inline
void
eye(ArmaMat<eT>& X, const ullint_t n)
{
    X = arma::eye<ArmaMat<eT>>(n,n);
}
#endif

#ifdef STATS_ENABLE_BLAZE_WRAPPERS
template<typename eT, bool To>
statslib_inline
void
eye(BlazeMat<eT,To>& X, const ullint_t n)
{
    X.resize(n,n);
    X = eT(0.0);
    for (ullint_t j=0; j < n; ++j) {
        X(j,j) = eT(1.0);
    }
}
#endif

#ifdef STATS_ENABLE_EIGEN_WRAPPERS
template<typename eT, int iTr, int iTc>
statslib_inline
void
eye(EigenMat<eT,iTr,iTc>& X, const ullint_t n)
{
    X.resize(n,n);
    X.setIdentity();
}
#endif
