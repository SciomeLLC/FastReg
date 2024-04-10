#' Randomly generated binary phenotype data used for running vignettes
#' A compressed text file with 500 individuals x 1 binary response
#' @format A data frame with 500 rows and 2 columns:
#' \describe{
#'   \item{ind}{\code{character} A unique alphanumeric identifier for each individual. Example: "IND123".}
#'   \item{bin.resp}{\code{integer} Binary phenotype data where 1 indicates the presence of a trait and 0 indicates its absence.}
#' }
"testdata_500_by_1k.bin.pheno"

#' Randomly generated covariate data used for running vignettes
#' A compressed text file with 500 individuals x 3 covariates (age, sex, treatment)
#' @format A data frame with 500 rows and 2 columns:
#' \describe{
#'   \item{ind}{\code{character} A unique alphanumeric identifier for each individual. Example: "IND123".}
#'   \item{age}{\code{integer} The age of the individual.}
#'   \item{sex}{\code{character} The sex of the individual - either M or F.}
#'   \item{treatment}{\code{character} The treatment that was administered to the individual. Test or Placebo.}
#' }
"testdata_500_by_1k.covar"


#' Randomly generated numerical phenotye data used for running vignettes
#' A compressed text file with 500 individuals x 1 numerical response
#' @format A data frame with 500 rows and 2 columns:
#' \describe{
#'   \item{ind}{\code{character} A unique alphanumeric identifier for each individual. Example: "IND123".}
#'   \item{num.resp}{\code{integer} Numerical phenotype data.}
#' }
"testdata_500_by_1k.num.pheno"

#' Randomly generated binary POI data used for running vignettes
#' A compressed text file with 500 individuals x 1000 binary responses
#' @format A data frame with 500 rows and 1000 columns:
#' \describe{
#'   \item{ind}{\code{character} A unique alphanumeric identifier for each individual. Example: "IND123".}
#'   \item{rs0001}{\code{integer} Binary response for the 1st poi. May contain no response at all.}
#'   \item{rs0002}{\code{integer} Binary response for the 2nd poi. May contain no response at all.}
#'   \item{rs####}{\code{integer} Binary response for the nth poi. May contain no response at all.}
#' }
"testdata_500_by_1k.poi"