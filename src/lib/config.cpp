#include <string>
#include <iostream>
#include <covariate.h>
#include <config.h>

using namespace arma;

template<class T>
T Config::get(const std::string &key)
{
    auto it = values.find(key);
    if (it == values.end())
    {
        Rcpp::Rcerr << "Key not found in config file: " + key << std::endl;
    }
    T result;
    std::istringstream stream(it->second);
    if (!(stream >> result))
    {
        Rcpp::Rcerr << "Failed to parse value for key: " + key << std::endl;
    }
    return result;
}

bool Config::has_key(const std::string &key)
{
    return values.find(key) != values.end();
}

std::string Config::get_value(std::string key, std::string def)
{
    if (values.find(key) == values.end())
    {
        return def;
    }
    return values[key];
}

void Config::trim(std::string &s)
{
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
    s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
}

void Config::parse(std::string file_path, const char delim, const char comment)
{
    std::ifstream conf_file(file_path);
    if (!conf_file.is_open())
    {
        Rcpp::Rcerr << "Error: unable to open configuration file " << file_path << std::endl;
        return;
    }
    std::string line;
    while (getline(conf_file, line))
    {
        size_t comment_start = line.find(comment);
        if (comment_start != std::string::npos)
        {
            line.erase(comment_start);
        }
        size_t delimiter_pos = line.find(delim);
        if (delimiter_pos == std::string::npos)
        {
            continue;
        }
        std::string key = line.substr(0, delimiter_pos);
        std::string value = line.substr(delimiter_pos + 1);
        trim(key);
        trim(value);
        if (key.empty() || value.empty())
        {
            continue;
        }
        values[key] = value;
    }
    conf_file.close();
}

void Config::validate_keys()
{
    std::vector<std::string> required_keys = {
        "pheno.file",
        "covar.file",
        "POI.file",
        "output.dir",
        "pheno.rowname.cols",
        "covar.rowname.cols",
        "phenotype",
        "regression.type"};

    for (const auto &key : required_keys)
    {
        if (!has_key(key))
        {
            Rcpp::Rcerr << "Error: configuration file is missing required key '" << key << "'" << std::endl;
        }
    }
}

void Config::validate_required_files()
{
    std::vector<std::string> file_keys = {
        "pheno.file",
        "covar.file",
        "POI.file"};

    for (const auto &key : file_keys)
    {
        std::string filename = values[key];
        if (!fs::exists(filename))
        {
            Rcpp::Rcerr << "Error: file '" << key << "' does not exist: " << filename << "\nPlease check that the path is correct." << std::endl;
        }
    }
}
void Config::set_default_values()
{
    std::unordered_map<std::string, std::string> defaults{
        {"no.intercept", "0"},
        {"maf.threshold", "0.01"},
        {"hwe.threshold", "1e-14"},
        {"colinearity.rsq", "0.99"},
        {"POI.file.format", "txt"},
        {"POI.file.delim", "tab"},
        {"pheno.file.delim", "tab"},
        {"covar.file.delim", "tab"},
        {"POI.type", "genotypes"},
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
        {"max.cores", "-1"}};

    for (auto item : defaults)
    {
        if (!has_key(item.first))
        {
            values[item.first] = item.second;
        }
    }

    if (!has_key("POI.effect.type"))
    {
        if (values.at("POI.type") == "genotype")
        {
            values["POI.effect.type"] = "additive";
        }
        else
        {
            values["POI.effect.type"] = "dosage";
        }
    }
}

void Config::validate_args()
{
    std::string empty_str = "";
    std::string temp = get<std::string>("POI.file.format");
    if (!(temp == "txt" || temp == "plink" || temp == "h5"))
    {
        throw std::invalid_argument("pheno.file.format not supported");
    }
    POI_file_format = temp;

    temp.assign(get<std::string>("POI.file.delim"));
    if (POI_file_format == "txt" && !(temp == "tab" || temp == "comma"))
    {
        throw std::invalid_argument("POI.file.delim not supported");
    }
    POI_file_delim = temp;

    temp.assign(get<std::string>("POI.effect.type"));
    if (!(temp == "dosage" || temp == "additive" || temp == "recessive" || temp == "dominant"))
    {
        throw std::invalid_argument("invalid POI.effect.type");
    }
    POI_effect_type = temp;

    temp.assign(get<std::string>("regression.type"));
    if (!(temp == "logistic" || temp == "linear"))
    {
        throw std::invalid_argument("invalid regression.type");
    }
    regression_type = temp;

    temp.assign(get<std::string>("pheno.file.delim"));
    if (!(temp == "tab" || temp == "comma" || temp == "semicolon"))
    {
        throw std::invalid_argument("pheno.file.delim not supported. Valid values are tab, comma and semicolon.");
    }
    pheno_file_delim = temp;

    temp.assign(get<std::string>("covar.file.delim"));
    if (!(temp == "tab" || temp == "comma" || temp == "semicolon"))
    {
        throw std::invalid_argument("cov.file.delim not supported. Valid values are tab, comma and semicolon.");
    }
    covar_file_delim = temp;

    no_intercept = get<bool>("no.intercept");
    output_exclude_covar = get<bool>("output.exclude.covar");
    verbose = get<bool>("verbose");
    compress_results = get<bool>("compress.results");

    POI_type = get<std::string>("POI.type");
    maf_threshold = get<double>("maf.threshold");
    if (maf_threshold > 0.5 || maf_threshold < 0)
    {
        throw std::invalid_argument("maf.threshold out of conventional bound");
    }

    hwe_threshold = get<double>("hwe.threshold");
    if (hwe_threshold > 0.5 || hwe_threshold < 0)
    {
        throw std::invalid_argument("hwe.threshold out of conventional bound");
    }

    colinearity_rsq = get<double>("colinearity.rsq");
    if (colinearity_rsq < 0.8 || colinearity_rsq > 1)
    {
        throw std::invalid_argument("colinearity.rsq out of conventional bound");
    }

    max_iter = get<int>("max.iter");
    max_cores = get<int>("max.cores");
    poi_block_size = get<int>("poi.block.size");

    temp.assign(get<std::string>("output.file.format"));
    if (!(temp == "long" || temp == "wide" || temp == "specific"))
    {
        throw std::invalid_argument("invalid output.file.format");
    }
    output_file_format = temp;

    temp.assign(get<std::string>("Pvalue.type"));
    if (!(temp == "t.dist" || temp == "norm.dist"))
    {
        throw std::invalid_argument("Pvalue.type must be t.dist or norm.dist");
    }
    p_value_type = temp;

    // required values
    pheno_file = values.at("pheno.file");
    covar_file = values.at("covar.file");
    POI_file = values.at("POI.file");
    output_dir = values.at("output.dir");
    phenotype = values.at("phenotype");
    pheno_rowname_cols = values.at("pheno.rowname.cols");
    covar_rowname_cols = values.at("covar.rowname.cols");

    split_by = split(values.at("split.by"), ",", "", 0);
    POI_covar_interactions = split(values.at("POI.covar.interactions"), ",", "", 0);
    covariate_terms = values.at("covariate.terms");

    // Subject subset
    subject_subset_rowname_cols = get_value("subject.subset.rowname.cols", empty_str);
    subject_subset_delim = get_value("subject.subset.delim", empty_str);
    subject_subset_file = get_value("subject.subset.file", empty_str);
    if (subject_subset_file != "")
    {
        if (!fs::exists(subject_subset_file))
        {
            Rcpp::Rcerr << "Error: file does not exist: " << subject_subset_file << "\nPlease check that the path is correct." << std::endl;
        }
    }

    // POI subset
    POI_subset_file = get_value("POI.subset.file", empty_str);
    POI_subset_file_delim = get_value("POI.subset.file.delim", empty_str);
    POI_subset_rowname_col = get_value("POI.subset.rowname.cols", empty_str);

    if (POI_subset_file != "")
    {
        if (!fs::exists(POI_subset_file))
        {
            Rcpp::Rcerr << "Error: file does not exist: " << POI_subset_file << "\nPlease check that the path is correct." << std::endl;
        }
    }

    // Optional covariate parameters
    covariates = split(get_value("covariates", empty_str), ",", "", 0);
    size_t cov_size = covariates.size();
    covariate_type = split(get_value("covariate.type", empty_str), ",", "", cov_size);
    covariate_standardize = split(get_value("covariate.standardize", empty_str), ",", "0", cov_size);
    covariate_levels = split(get_value("covariate.levels", empty_str), ";", "", cov_size);
    covariate_ref_level = split(get_value("covariate.ref.level", empty_str), ",", "", cov_size);


    if (cov_size == 1 && covariates[0].empty()) {
        cov_size = 0;
    }
    if (cov_size > 0) {
        if (!std::all_of(
                covariate_type.begin(),
                covariate_type.end(),
                [](const std::string& kv) { return kv == "numeric" || kv == "categorical" || kv == "count"; }
            )
        ) {
            Rcpp::stop("invalid covariate.type");
        }

        std::vector<int> covariate_standardize_ints;

        for(auto& it : covariate_standardize) {
            covariate_standardize_ints.push_back(std::stoi(it));
        }

        
        if(covariate_type.size() != cov_size) {
            Rcpp::stop("number of covariate.type != number of covariates.");
        }

        if(covariate_standardize.size() != cov_size) {
            Rcpp::stop("number of covariate.standardize != number of covariates.");
        }

        if(covariate_ref_level.size() != cov_size) {
            Rcpp::stop("number of covariate.ref.level != number of covariates.");
        }

        for(unsigned int i = 0; i < cov_size; i++) {
            Covariate c(covariates[i], covariate_type[i], covariate_ref_level[i], covariate_levels[i], covariate_standardize_ints[i]);
            covs.push_back(c);
        }
    }
}

std::vector<std::string> Config::split(std::string val, std::string delim, std::string default_str_val, unsigned int size) {
    std::vector<std::string> split_result;

    size_t pos = 0;
    std::string token;

    int count = 0;
    while((pos = val.find(delim)) != std::string::npos) {
        token = val.substr(0, pos);
        split_result.push_back(token);
        val.erase(0, pos + delim.length());
        count++;
    }
    split_result.push_back(val);
    
    while(split_result.size() < size) {
        split_result.push_back(default_str_val);
    }

    return split_result;
}
