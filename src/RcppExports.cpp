// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// FastRegImportCpp
int FastRegImportCpp(std::string dataFile, std::string h5File, int headerRow, int idCol, int dataCol, float buffSize, bool transpose, int chunkEdge, bool vcf, std::string delim, bool gz, int poiPerFile, bool singleFile, int serverThreads, float serverMem);
RcppExport SEXP _FastReg_FastRegImportCpp(SEXP dataFileSEXP, SEXP h5FileSEXP, SEXP headerRowSEXP, SEXP idColSEXP, SEXP dataColSEXP, SEXP buffSizeSEXP, SEXP transposeSEXP, SEXP chunkEdgeSEXP, SEXP vcfSEXP, SEXP delimSEXP, SEXP gzSEXP, SEXP poiPerFileSEXP, SEXP singleFileSEXP, SEXP serverThreadsSEXP, SEXP serverMemSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type dataFile(dataFileSEXP);
    Rcpp::traits::input_parameter< std::string >::type h5File(h5FileSEXP);
    Rcpp::traits::input_parameter< int >::type headerRow(headerRowSEXP);
    Rcpp::traits::input_parameter< int >::type idCol(idColSEXP);
    Rcpp::traits::input_parameter< int >::type dataCol(dataColSEXP);
    Rcpp::traits::input_parameter< float >::type buffSize(buffSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type transpose(transposeSEXP);
    Rcpp::traits::input_parameter< int >::type chunkEdge(chunkEdgeSEXP);
    Rcpp::traits::input_parameter< bool >::type vcf(vcfSEXP);
    Rcpp::traits::input_parameter< std::string >::type delim(delimSEXP);
    Rcpp::traits::input_parameter< bool >::type gz(gzSEXP);
    Rcpp::traits::input_parameter< int >::type poiPerFile(poiPerFileSEXP);
    Rcpp::traits::input_parameter< bool >::type singleFile(singleFileSEXP);
    Rcpp::traits::input_parameter< int >::type serverThreads(serverThreadsSEXP);
    Rcpp::traits::input_parameter< float >::type serverMem(serverMemSEXP);
    rcpp_result_gen = Rcpp::wrap(FastRegImportCpp(dataFile, h5File, headerRow, idCol, dataCol, buffSize, transpose, chunkEdge, vcf, delim, gz, poiPerFile, singleFile, serverThreads, serverMem));
    return rcpp_result_gen;
END_RCPP
}
// FastRegCpp
void FastRegCpp(const std::string phenotype, const std::string regression_type, const std::string pvalue_dist, bool output_exclude_covar, double maf_threshold, double hwe_threshold, bool no_intercept, double colinearity_rsq, int poi_block_size, int max_iter, double rel_conv_tolerance, double abs_conv_tolderance, int max_openmp_threads, const std::string pheno_file, const std::string pheno_rowname_cols, const std::string pheno_file_delim, const std::string covar_file, const std::string covar_rowname_cols, const std::string covar_file_delim, const std::string poi_file_dir, const std::string poi_file_delim, const std::string poi_file_format, const std::string poi_type, const std::string poi_effect_type, const Rcpp::StringVector covariates, const Rcpp::StringVector covariate_type, const Rcpp::LogicalVector covariate_standardize, const Rcpp::StringVector covariate_levels, const Rcpp::StringVector covariate_ref_level, const Rcpp::StringVector POI_covar_interactions_str, const Rcpp::StringVector split_by_str, const std::string output_dir, bool compress_results, int max_workers);
RcppExport SEXP _FastReg_FastRegCpp(SEXP phenotypeSEXP, SEXP regression_typeSEXP, SEXP pvalue_distSEXP, SEXP output_exclude_covarSEXP, SEXP maf_thresholdSEXP, SEXP hwe_thresholdSEXP, SEXP no_interceptSEXP, SEXP colinearity_rsqSEXP, SEXP poi_block_sizeSEXP, SEXP max_iterSEXP, SEXP rel_conv_toleranceSEXP, SEXP abs_conv_tolderanceSEXP, SEXP max_openmp_threadsSEXP, SEXP pheno_fileSEXP, SEXP pheno_rowname_colsSEXP, SEXP pheno_file_delimSEXP, SEXP covar_fileSEXP, SEXP covar_rowname_colsSEXP, SEXP covar_file_delimSEXP, SEXP poi_file_dirSEXP, SEXP poi_file_delimSEXP, SEXP poi_file_formatSEXP, SEXP poi_typeSEXP, SEXP poi_effect_typeSEXP, SEXP covariatesSEXP, SEXP covariate_typeSEXP, SEXP covariate_standardizeSEXP, SEXP covariate_levelsSEXP, SEXP covariate_ref_levelSEXP, SEXP POI_covar_interactions_strSEXP, SEXP split_by_strSEXP, SEXP output_dirSEXP, SEXP compress_resultsSEXP, SEXP max_workersSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type regression_type(regression_typeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type pvalue_dist(pvalue_distSEXP);
    Rcpp::traits::input_parameter< bool >::type output_exclude_covar(output_exclude_covarSEXP);
    Rcpp::traits::input_parameter< double >::type maf_threshold(maf_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type hwe_threshold(hwe_thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type no_intercept(no_interceptSEXP);
    Rcpp::traits::input_parameter< double >::type colinearity_rsq(colinearity_rsqSEXP);
    Rcpp::traits::input_parameter< int >::type poi_block_size(poi_block_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type rel_conv_tolerance(rel_conv_toleranceSEXP);
    Rcpp::traits::input_parameter< double >::type abs_conv_tolderance(abs_conv_tolderanceSEXP);
    Rcpp::traits::input_parameter< int >::type max_openmp_threads(max_openmp_threadsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type pheno_file(pheno_fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type pheno_rowname_cols(pheno_rowname_colsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type pheno_file_delim(pheno_file_delimSEXP);
    Rcpp::traits::input_parameter< const std::string >::type covar_file(covar_fileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type covar_rowname_cols(covar_rowname_colsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type covar_file_delim(covar_file_delimSEXP);
    Rcpp::traits::input_parameter< const std::string >::type poi_file_dir(poi_file_dirSEXP);
    Rcpp::traits::input_parameter< const std::string >::type poi_file_delim(poi_file_delimSEXP);
    Rcpp::traits::input_parameter< const std::string >::type poi_file_format(poi_file_formatSEXP);
    Rcpp::traits::input_parameter< const std::string >::type poi_type(poi_typeSEXP);
    Rcpp::traits::input_parameter< const std::string >::type poi_effect_type(poi_effect_typeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type covariate_type(covariate_typeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector >::type covariate_standardize(covariate_standardizeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type covariate_levels(covariate_levelsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type covariate_ref_level(covariate_ref_levelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type POI_covar_interactions_str(POI_covar_interactions_strSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type split_by_str(split_by_strSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output_dir(output_dirSEXP);
    Rcpp::traits::input_parameter< bool >::type compress_results(compress_resultsSEXP);
    Rcpp::traits::input_parameter< int >::type max_workers(max_workersSEXP);
    FastRegCpp(phenotype, regression_type, pvalue_dist, output_exclude_covar, maf_threshold, hwe_threshold, no_intercept, colinearity_rsq, poi_block_size, max_iter, rel_conv_tolerance, abs_conv_tolderance, max_openmp_threads, pheno_file, pheno_rowname_cols, pheno_file_delim, covar_file, covar_rowname_cols, covar_file_delim, poi_file_dir, poi_file_delim, poi_file_format, poi_type, poi_effect_type, covariates, covariate_type, covariate_standardize, covariate_levels, covariate_ref_level, POI_covar_interactions_str, split_by_str, output_dir, compress_results, max_workers);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastReg_FastRegImportCpp", (DL_FUNC) &_FastReg_FastRegImportCpp, 15},
    {"_FastReg_FastRegCpp", (DL_FUNC) &_FastReg_FastRegCpp, 34},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
