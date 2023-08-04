// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <chrono>
#define R_NO_REMAP
#include <Rmath.h>
#include <fr_matrix.h>
#include <covariate.h>
#include <stats.hpp>
#include <regression.h>

#undef pnorm
using namespace arma;

arma::colvec RegressionBase::t_dist(arma::colvec abs_z, int df, bool log_form) {
    return -1*(stats2::pt(abs_z, df, log_form) + log(2))/log(10);
}

arma::colvec RegressionBase::norm_dist(arma::colvec abs_z, int df, bool log_form) {
    return -1*(arma::log_normpdf(abs_z) + log(2))/log(10);
}

void LogisticRegression::run(
    FRMatrix& cov, 
    FRMatrix& pheno, 
    FRMatrix& poi_data,
    FRMatrix& interactions,
    arma::umat& W2,
    FRMatrix& beta_est,
    FRMatrix& se_beta,
    FRMatrix& neglog10_pvl,
    int max_iter, 
    bool is_t_dist) {
    arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;

    // create a pointer to the specified distribution function
    arma::colvec (*dist_func)(arma::colvec, int, bool) = is_t_dist == true ? t_dist : norm_dist;
    // arma::colvec (*dist_func)(arma::colvec, int) = t_dist;
    
    #pragma omp parallel for
    for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++) {
        arma::mat A(n_parms, n_parms, arma::fill::zeros);
        // arma::mat poi_mat = arma::conv_to<arma::mat>::from(poi_data.data.col(poi_col));

        // Initialize beta
        arma::colvec beta(n_parms, arma::fill::zeros);
        arma::mat cov_w_mat = cov.data;
        arma::mat int_w_mat = interactions.data; 
        arma::ucolvec w2_col = W2.col(poi_col);
        // int_w_mat.each_col([&poi_data, &w2_col, &poi_col](arma::vec& a){(a%poi_data.data.col(poi_col))%w2_col;});
        // int_w_mat.each_col() %= poi_data.data.col(poi_col) % w2_col;
        int_w_mat = interactions.data % arma::repmat(poi_data.data.col(poi_col), 1, interactions.data.n_cols);
        arma::uword first_chunk = cov_w_mat.n_cols - 1;
        arma::uword second_chunk = n_parms - 1;
        arma::span first = arma::span(0, first_chunk);
        arma::span second = arma::span(first_chunk + 1, second_chunk);
        // arma::span col_1 = arma::span(0,0);

        
        for (int iter = 0; iter < max_iter; iter++) {
            arma::mat eta(cov.data.n_rows, 1, arma::fill::zeros);
            eta += cov_w_mat * beta.subvec(first);
            eta += int_w_mat * beta.subvec(second);
            
            arma::mat p = 1 / (1 + arma::exp(-eta));
            arma::mat W1 = p % (1-p) % w2_col;
            
            // cov_w_mat.each_col([&W1](arma::vec &a){a%W1;});
            // int_w_mat.each_col([&W1](arma::vec &a){a%W1;});
            // cov_w_mat.each_col() %= W1;
            // int_w_mat.each_col() %= W1;
            arma::mat temp1 = cov_w_mat % arma::repmat(W1, 1, cov_w_mat.n_cols);
            arma::mat temp2 = int_w_mat % arma::repmat(W1, 1, int_w_mat.n_cols);


            A.submat(first, first) = temp1.t() * cov_w_mat; // C'C
            A.submat(second, second) = temp2.t() * int_w_mat; // I'I
            A.submat(first, second) = temp1.t() * int_w_mat; // C'I
            A.submat(second, first) = A.submat(first, second).t(); // I'C

            arma::mat z = w2_col % (pheno.data - p);

            arma::mat B = arma::mat(n_parms, 1, arma::fill::zeros);
            B.submat(first, arma::span(0, 0)) = cov_w_mat.t()*z;
            B.submat(second, arma::span(0, 0)) = int_w_mat.t()*z;

            beta = beta + arma::solve(A, B);
        }
        
        int df = arma::as_scalar(arma::sum(w2_col, 0)) - n_parms;
        arma::colvec temp_se = arma::sqrt(arma::abs(arma::diagvec(arma::pinv(A))));

        beta_est.data.col(poi_col) = beta;
        se_beta.data.col(poi_col) = arma::sqrt(arma::abs(arma::diagvec(arma::pinv(A))));
        neglog10_pvl.data.col(poi_col) = (*dist_func)(arma::abs(beta/temp_se), df, true);
    }
}

void LinearRegression::run(
    FRMatrix& cov, 
    FRMatrix& pheno, 
    FRMatrix& poi_data,
    FRMatrix& interactions,
    arma::umat& W2,
    FRMatrix& beta_est,
    FRMatrix& se_beta,
    FRMatrix& neglog10_pvl,
    int max_iter, 
    bool is_t_dist) {

    arma::uword n_parms = cov.data.n_cols + interactions.data.n_cols;
    arma::span col_1 = arma::span(0,0);
    // arma::colvec (*dist_func)(Rcpp::NumericVector, int) = is_t_dist == true ? &t_dist : &norm_dist;
    
    arma::colvec (*dist_func)(arma::colvec, int, bool) = is_t_dist == true ? t_dist : norm_dist;
    #pragma omp parallel for
    for (arma::uword poi_col = 0; poi_col < poi_data.data.n_cols; poi_col++) {
        // Initialize beta
        arma::mat beta(n_parms, 1, arma::fill::zeros);
        arma::mat cov_w_mat = cov.data;
        arma::colvec poi_mat = poi_data.data.col(poi_col);
        arma::ucolvec w2_col = W2.col(poi_col);
        // arma::mat int_w_mat = interactions.data;
        // int_w_mat.each_col([&poi_mat, &w2_col](arma::vec& a){(a%poi_mat)%w2_col;});
        arma::mat int_w_mat = interactions.data % arma::repmat(poi_mat, 1, interactions.data.n_cols);
        arma::uword first_chunk = cov_w_mat.n_cols - 1;
        arma::uword second_chunk = n_parms - 1;
        arma::span first = arma::span(0, first_chunk);
        arma::span second = arma::span(first_chunk + 1, second_chunk);
        // cov_w_mat.each_col([&w2_col](arma::vec &a){a%w2_col;});
        // int_w_mat.each_col([&w2_col](arma::vec &a) {a%w2_col;});
        arma::mat temp1 = cov_w_mat % arma::repmat(w2_col, 1, cov_w_mat.n_cols);
        arma::mat temp2 = int_w_mat % arma::repmat(w2_col, 1, int_w_mat.n_cols);
        
        arma::mat A(n_parms, n_parms, arma::fill::zeros);
        A.submat(first, first) = temp1.t() * cov_w_mat; // C'C
        A.submat(second, second) = temp2.t() * int_w_mat; // I'I
        A.submat(first, second) = temp1.t() * int_w_mat; // C'I
        A.submat(second, first) = A.submat(first, second).t(); // I'C

        arma::mat z = w2_col % (pheno.data);

        arma::mat B = arma::mat(n_parms, 1, arma::fill::zeros);
        B.submat(first, col_1) = cov_w_mat.t()*z;
        B.submat(second, col_1) = int_w_mat.t()*z;
        beta = arma::solve(A, B);

        beta_est.data.col(poi_col) = beta;
        arma::mat eta(cov.data.n_rows, 1, arma::fill::zeros);
        eta += cov_w_mat * beta.submat(first, col_1);
        eta += int_w_mat * beta.submat(second, col_1);
        int df = arma::as_scalar(arma::sum(W2.col(poi_col), 0)) - n_parms;
        double mse = arma::conv_to<double>::from(W2.col(poi_col).t() * arma::square(pheno.data - eta))/(double)df;

        se_beta.data.col(poi_col) = arma::sqrt(mse * arma::abs(arma::diagvec(arma::pinv(A))));
        arma::mat temp_se = se_beta.data.col(poi_col);
        neglog10_pvl.data.col(poi_col) = dist_func(arma::abs(beta/temp_se), df, true);
    }
}