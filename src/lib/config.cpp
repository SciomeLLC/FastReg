#include <config.h>
#include <covariate.h>
#include <iostream>
#include <string>

using namespace arma;

template <class T> T Config::get(const std::string &key) {
  auto it = values.find(key);
  if (it == values.end()) {
    Rcpp::Rcerr << "Key not found in config file: " + key << std::endl;
  }
  T result;
  std::istringstream stream(it->second);
  if (!(stream >> result)) {
    Rcpp::Rcerr << "Failed to parse value for key: " + key << std::endl;
  }
  return result;
}

bool Config::has_key(const std::string &key) {
  return values.find(key) != values.end();
}

std::string Config::get_value(std::string key, std::string def) {
  if (values.find(key) == values.end()) {
    return def;
  }
  return values[key];
}

void Config::trim(std::string &s) {
  s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
  s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
}

void Config::parse(std::string file_path, const char delim,
                   const char comment) {
  std::ifstream conf_file(file_path);
  if (!conf_file.is_open()) {
    Rcpp::Rcerr << "Error: unable to open configuration file " << file_path
                << std::endl;
    return;
  }
  std::string line;
  while (getline(conf_file, line)) {
    size_t comment_start = line.find(comment);
    if (comment_start != std::string::npos) {
      line.erase(comment_start);
    }
    size_t delimiter_pos = line.find(delim);
    if (delimiter_pos == std::string::npos) {
      continue;
    }
    std::string key = line.substr(0, delimiter_pos);
    std::string value = line.substr(delimiter_pos + 1);
    trim(key);
    trim(value);
    if (key.empty() || value.empty()) {
      continue;
    }
    values[key] = value;
  }
  conf_file.close();
}

void Config::validate_keys() {
  std::vector<std::string> required_keys = {
      "pheno.file", "covar.file",         "POI.file",
      "output.dir", "pheno.rowname.cols", "covar.rowname.cols",
      "phenotype",  "regression.type"};

  for (const auto &key : required_keys) {
    if (!has_key(key)) {
      Rcpp::Rcerr << "Error: configuration file is missing required key '"
                  << key << "'" << std::endl;
    }
  }
}

std::vector<std::string>
Config::convert_stringV_to_string_arr(Rcpp::StringVector stringV) {
  std::vector<std::string> vstrings(stringV.size());
  for (int i = 0; i < stringV.size(); i++) {
    vstrings[i] = stringV(i);
  }

  return vstrings;
}
void Config::validate_covariate_config(
    Rcpp::StringVector covariates_str, Rcpp::StringVector covariate_type_str,
    Rcpp::LogicalVector covariate_standardize_str,
    Rcpp::StringVector covariate_levels_str,
    Rcpp::StringVector covariate_ref_level_str,
    Rcpp::StringVector POI_covar_interactions_str) {
  // Optional covariate parameters
  // covariates = split(covariates_str, ",", "", 0);
  covariates = convert_stringV_to_string_arr(covariates_str);
  size_t cov_size = covariates.size();
  covariate_type = convert_stringV_to_string_arr(covariate_type_str);

  // covariate_standardize = split(covariate_standardize_str, ",", "0",
  // cov_size);
  covariate_levels = convert_stringV_to_string_arr(covariate_levels_str);
  covariate_ref_level = convert_stringV_to_string_arr(covariate_ref_level_str);

  for (const int &cov_std : covariate_standardize_str) {
    covariate_standardize.push_back((bool)cov_std);
  }

  // split_by = split(split_by_str, ",", "", 0);
  POI_covar_interactions =
      convert_stringV_to_string_arr(POI_covar_interactions_str);
  if (cov_size == 1 && covariates[0].empty()) {
    cov_size = 0;
  }
  if (cov_size > 0) {
    if (!std::all_of(covariate_type.begin(), covariate_type.end(),
                     [](const std::string &kv) {
                       return kv == "numeric" || kv == "categorical" ||
                              kv == "count";
                     })) {
      Rcpp::stop("invalid covariate.type");
    }

    if (covariate_type.size() != cov_size) {
      Rcpp::stop("number of covariate.type != number of covariates.");
    }

    if (covariate_standardize.size() != cov_size) {
      Rcpp::stop("number of covariate.standardize != number of covariates.");
    }

    if (covariate_levels.size() != cov_size) {
      Rcpp::stop("number of covariate.levels != number of covariates.");
    }

    if (covariate_ref_level.size() != cov_size) {
      Rcpp::stop("number of covariate.ref.level != number of covariates.");
    }

    for (unsigned int i = 0; i < cov_size; i++) {
      Covariate c(covariates[i], covariate_type[i], covariate_ref_level[i],
                  covariate_levels[i], covariate_standardize_str[i]);
      covs.push_back(c);
    }
  }
}
void Config::validate_required_files() {
  // std::vector<std::string> file_keys = {
  //     pheno_file,
  //     covar_file,
  //     POI_file};

  if (!fs::exists(pheno_file)) {
    Rcpp::Rcerr << "Error: file does not exist: " << pheno_file
                << "\nPlease check that the path is correct." << std::endl;
  }
  if (!fs::exists(covar_file)) {
    Rcpp::Rcerr << "Error: file does not exist: " << covar_file
                << "\nPlease check that the path is correct." << std::endl;
  }
  if (!fs::exists(POI_file_dir)) {
    Rcpp::Rcerr << "Error: Directory does not exist - " << POI_file_dir
                << "\nPlease check that the path is correct." << std::endl;
  }
}

void Config::set_default_values() {
  std::unordered_map<std::string, std::string> defaults{
      {"no.intercept", "0"},
      {"maf.threshold", "1e-13"},
      {"hwe.threshold", "1e-13"},
      {"colinearity.rsq", "1.0"},
      {"POI.file.format", "txt"},
      {"POI.file.delim", "tab"},
      {"pheno.file.delim", "tab"},
      {"covar.file.delim", "tab"},
      {"POI.type", "genotype"},
      {"POI.subset.file", ""},
      {"POI.covar.interactions", ""},
      {"split.by", ""},
      {"output.file.format", "long"},
      {"output.exclude.covar", "0"},
      {"poi.block.size", "0"},
      {"covariate.terms", ""},
      {"max.iter", "6"},
      {"Pvalue.type", "t.dist"},
      {"verbose", "1"},
      {"compress.results", "1"},
      {"max.threads", "-1"},
      {"rel.conv.tolerance", "0.0001"},
      {"abs.conv.tolerance", "0.0001"}};

  for (auto item : defaults) {
    if (!has_key(item.first)) {
      values[item.first] = item.second;
    }
  }

  if (!has_key("POI.effect.type")) {
    if (values.at("POI.type") == "genotype") {
      values["POI.effect.type"] = "additive";
    } else {
      values["POI.effect.type"] = "dosage";
    }
  }
}

void Config::validate_args() {
  std::string empty_str = "";
  if (!(POI_file_format == "txt" || POI_file_format == "plink" ||
        POI_file_format == "h5" || POI_file_format == "bed")) {
    Rcpp::stop("POI.file.format not supported");
  }

  if (POI_file_format == "txt" &&
      !(POI_file_delim == "tab" || POI_file_delim == "comma")) {
    Rcpp::stop("POI.file.delim not supported");
  }

  if (!(POI_effect_type == "dosage" || POI_effect_type == "additive" ||
        POI_effect_type == "recessive" || POI_effect_type == "dominant")) {
    Rcpp::stop("invalid POI.effect.type");
  }

  if (!(regression_type == "logistic" || regression_type == "linear")) {
    Rcpp::stop("invalid regression.type");
  }

  if (!(pheno_file_delim == "tab" || pheno_file_delim == "comma" ||
        pheno_file_delim == "semicolon")) {
    Rcpp::stop("pheno.file.delim not supported. Valid values "
               "are tab, comma and semicolon.");
  }

  if (!(covar_file_delim == "tab" || covar_file_delim == "comma" ||
        covar_file_delim == "semicolon")) {
    Rcpp::stop("cov.file.delim not supported. Valid values "
               "are tab, comma and semicolon.");
  }

  if (maf_threshold > 0.5 || maf_threshold < 0) {
    Rcpp::stop("maf.threshold out of conventional bound");
  }

  if (hwe_threshold > 0.5 || hwe_threshold < 0) {
    Rcpp::stop("hwe.threshold out of conventional bound");
  }

  if (colinearity_rsq < 0.8 || colinearity_rsq > 1) {
    Rcpp::stop("colinearity.rsq out of conventional bound");
  }

  if (!(p_value_type == "t.dist" || p_value_type == "norm.dist")) {
    Rcpp::stop("Pvalue.type must be t.dist or norm.dist");
  }
}

std::vector<std::string> Config::split(std::string val, std::string delim,
                                       std::string default_str_val,
                                       unsigned int size) {
  std::vector<std::string> split_result;

  size_t pos = 0;
  std::string token;

  int count = 0;
  while ((pos = val.find(delim)) != std::string::npos) {
    token = val.substr(0, pos);
    split_result.push_back(token);
    val.erase(0, pos + delim.length());
    count++;
  }
  split_result.push_back(val);

  while (split_result.size() < size) {
    split_result.push_back(default_str_val);
  }

  return split_result;
}
// get a list of poi files from the poi_file_dir directory with the
// poi_file_type extension
void Config::get_poi_files() {
  std::string extension;
  if (POI_file_format == "bed") {
    extension = ".bed";
  } else {
    extension = ".h5";
  }

  for (auto &entry : fs::directory_iterator(POI_file_dir)) {
    if (entry.path().extension() == extension) {
      poi_files.push_back(fs::absolute(entry.path()).string());
    }
  }
}

void Config::print() {
  Rcpp::Rcout << "-----------------------------------------" << std::endl;
  Rcpp::Rcout << "Running FastReg with configuration: " << std::endl;
  Rcpp::Rcout << "phenotype: " << phenotype << std::endl;
  Rcpp::Rcout << "regression_type: " << regression_type << std::endl;
  Rcpp::Rcout << "pvalue_dist: " << p_value_type << std::endl;
  Rcpp::Rcout << "output_exclude_covar: " << output_exclude_covar << std::endl;
  Rcpp::Rcout << "maf_threshold: " << maf_threshold << std::endl;
  Rcpp::Rcout << "hwe_threshold: " << hwe_threshold << std::endl;
  Rcpp::Rcout << "no_intercept: " << no_intercept << std::endl;
  Rcpp::Rcout << "colinearity_rsq: " << colinearity_rsq << std::endl;
  Rcpp::Rcout << "poi_block_size: " << poi_block_size << std::endl;
  Rcpp::Rcout << "max_iter: " << max_iter << std::endl;
  Rcpp::Rcout << "rel_conv_tolerance: " << rel_conv_tolerance << std::endl;
  Rcpp::Rcout << "abs_conv_tolerance: " << abs_conv_tolerance << std::endl;
  Rcpp::Rcout << "max_openmp_threads: " << max_openmp_threads << std::endl;
  Rcpp::Rcout << "pheno_file: " << pheno_file << std::endl;
  Rcpp::Rcout << "pheno_rowname_cols: " << pheno_rowname_cols << std::endl;
  Rcpp::Rcout << "pheno_file_delim: " << pheno_file_delim << std::endl;
  Rcpp::Rcout << "covar_file: " << covar_file << std::endl;
  Rcpp::Rcout << "covar_rowname_cols: " << covar_rowname_cols << std::endl;
  Rcpp::Rcout << "covar_file_delim: " << covar_file_delim << std::endl;
  Rcpp::Rcout << "poi_file_dir: " << POI_file_dir << std::endl;
  Rcpp::Rcout << "poi_file_delim: " << POI_file_delim << std::endl;
  Rcpp::Rcout << "poi_file_format: " << POI_file_format << std::endl;
  Rcpp::Rcout << "poi_type: " << POI_type << std::endl;
  Rcpp::Rcout << "poi_effect_type: " << POI_effect_type << std::endl;
  for (size_t i = 0; i < covs.size(); i++) {
    covs[i].print();
  }
  Rcpp::Rcout << "POI_covar_interactions_str: " << POI_covar_interactions_str
              << std::endl;
  Rcpp::Rcout << "split_by_str: " << split_by_str << std::endl;
  Rcpp::Rcout << "output_dir: " << output_dir << std::endl;
  Rcpp::Rcout << "compress_results: " << compress_results << std::endl;
  Rcpp::Rcout << "max_workers: " << max_workers << std::endl;
  Rcpp::Rcout << "-----------------------------------------" << std::endl;
}