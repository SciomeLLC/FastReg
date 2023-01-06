parseConfigFile <- function(file, delim="\t", comment.char="#", protected.args=NULL) {

  if(!is.null(protected.args)) {
      args <- protected.args;
      protected.args <- names(args);
  } else {
    args <- list();
    protected.args <- c();
  }

  txt <- readLines(file);
  txt <- txt[(txt!="") & regexpr(paste0("^", comment.char),txt)==-1 & regexpr(paste0("^", delim),txt)==-1]; ## Ignore commented and blank lines

  config.vals <- data.frame(do.call(rbind, strsplit(txt,split=delim)),stringsAsFactors=FALSE);
  colnames(config.vals) <- c("Argument", "Value");
  config.vals[,"Value"] <- gsub(" +$", "", gsub("^ +",   "", sub(pattern = "^\"", replacement = "", x = sub(pattern = "\"$", replacement = "", x = config.vals[,"Value"])))); ## trim leading and trailing blanks

  config.vals <- config.vals[!duplicated(config.vals[,"Argument"]),];
  rownames(config.vals) <- config.vals[,"Argument"];

  for(i in 1:nrow(config.vals))  {
    ival <- config.vals[i,"Value"];
    if(nchar(sub("^[.]", "", sub("^[-]", "", gsub("[0-9]","",ival))))==0) ival <- as.numeric(ival); #convert to number as necessary
    iname <- gsub(" ",".", config.vals[i,"Argument"]); #create valid varnames
    if(iname %in% c("covariates" , "covariate.type", "covariate.ref.level", "POI.Covar.interactions", "covariate.standardize", "Split.by")) {
      ival <- unlist(strsplit(ival, split=","));
    }

    if(iname == "covariate.standardize") {
      ival <- as.integer(ival);
      ival <- (ival==1)
    }

    if(iname == "covariate.levels") {
      ival <- unlist(strsplit(ival, split=";"));
      ival <- lapply(ival, FUN=function(x) unlist(strsplit(x, split=",")))
    }

    if(iname %in% c("POI.file.delim","pheno.file.delim", "covar.file.delim", "POI.subset.file.delim", "subject.subset.file.delim")) {
      if(ival == "tab") {
        ival <- "\t";
      } else if(ival=="comma") {
        ival <- ",";
      } else if(ival == "semicolon") {
        ival <- ";";
      }
    }

    na.ival <- (ival=="NA")
    if(any(na.ival)) ival[na.ival] <- NA;


    if(!(iname %in% protected.args)) args[[iname]] <- ival;
  }

  invisible(args);
}

assign.default.values <- function(args) {
  if(!is.null(args)) {
    protected.args <- names(args);
  } else {
    args <- list();
    protected.args <- c();
  }
  if(!("no.intercept" %in% protected.args)) args[["no.intercept"]] <- 0;
  if(!("maf.threshold" %in% protected.args)) args[["maf.threshold"]] <- 0.01;
  if(!("hwe.threshold" %in% protected.args)) args[["hwe.threshold"]] <- 0.05;
  if(!("colinearity.rsq" %in% protected.args)) args[["colinearity.rsq"]] <- 0.99;
  if(!("POI.file.format" %in% protected.args)) args[["POI.file.format"]] <- "txt";
  if(!("POI.file.delim" %in% protected.args)) args[["POI.file.delim"]] <- "\t";
  if(!("pheno.file.delim" %in% protected.args)) args[["pheno.file.delim"]] <-"\t";
  if(!("covar.file.delim" %in% protected.args)) args[["covar.file.delim"]] <- "\t";
  if(!("POI.type" %in% protected.args)) args[["POI.type"]]  <-"genotypes";
  if(!("subject.subset.file" %in% protected.args)) args[["subject.subset.file"]] <- NULL;
  if(!("POI.subset.file" %in% protected.args)) args[["POI.subset.file"]] <- NULL;
  if(!("POI.Covar.Interactions" %in% protected.args)) args[["POI.Covar.Interactions"]] <- NULL;
  if(!("Split.by" %in% protected.args)) args[["Split.by"]] <- NULL;
  if(!("output.dir" %in% protected.args)) {
      args[["output.dir"]] <- paste0("./",gsub("[.]txt", "",basename(args[["pheno.file"]])));
      dir.create(args[["output.dir"]], showWarnings=FALSE, recursive=TRUE);
  }
  if(!("outputfile.format" %in% protected.args)) args[["outputfile.format"]] <-"long";
  if(!("output.exclude.covar" %in% protected.args)) args[["output.exclude.covar"]] <- 0;
  if(!("poi.block.size" %in% protected.args)) args[["poi.block.size"]] <- 0;
  if(!("num.cores" %in% protected.args)) args[["num.cores"]] <- NULL;

  if(!("covariate.terms" %in% protected.args)) args[["covariate.terms"]] <- NULL;
  if(!("max.iter" %in% protected.args)) args[["max.iter"]] <- 6;


  if(!("POI.effect.type" %in% protected.args))   args[["POI.effect.type"]] <- ifelse(args[["POI.type"]]=="genotypes", "additive", "dosage");
  if(!("Pvalue.type" %in% protected.args)) args[["Pvalue.type"]] <- "t.dist";
  if(!("verbose" %in% protected.args)) args[["verbose"]] <- 0;
  invisible(args);
}

validate.args <- function(args) {
  required.args <- c("pheno.file", "covar.file", "POI.file", "output.dir", "pheno.rowname.cols","covar.rowname.cols",  "phenotype", "regression.type");
  missing.args <- required.args[!(required.args %in% names(args))];
  if(length(missing.args)>0) {
    stop("missing critical args:", paste(missing.args, collapse=", "));
  }

  if(!file.exists(args[["pheno.file"]])) stop("pheno.file does not exist");
  if(!file.exists(args[["covar.file"]])) stop("covar.file does not exist");
  if(!file.exists(args[["POI.file"]])) stop("POI.file does not exist");
  if(!is.null(args[["POI.subset.file"]])) {
    if(!file.exists(args[["POI.subset.file"]])) stop("POI.subset.file does not exist")
  }

  if(!is.null(args[["subject.subset.file"]])) {
    if(!file.exists(args[["subject.subset.file"]])) stop("subject.subset.file does not exist")
  }

  if(!dir.exists(args[["output.dir"]])) {
      if(!dir.create(args[["output.dir"]], recursive=TRUE, showWarnings=FALSE)) stop("output.dir does  could not be created")
  }

  if(!(args[["POI.file.format"]] %in% c("txt", "plink", "h5"))) stop("pheno.file.format not supported");
  if(args[["POI.file.format"]]=="txt" & !(args[["POI.file.delim"]] %in% c("\t", ","))) stop("POI.file.delim not supported")
  if(!(args[["POI.effect.type"]] %in% c("dosage", "additive", "recessive", "dominant"))) stop("invalid POI.effect.type")
  if(!(args[["regression.type"]] %in% c("logistic","linear"))) stop("invalid regression.type");

  if(!(args[["pheno.file.delim"]] %in% c("\t", ","))) stop("pheno.file.delim not supported")
  if(!(args[["covar.file.delim"]] %in% c("\t", ","))) stop("cov.file.delim not supported")

  if(!(args[["no.intercept"]] %in% c(0,1))) stop("no.intercept must be 0 or 1");
  if(args[["maf.threshold"]]>0.5 | args[["maf.threshold"]] < 0) stop("maf.threshold out off conventional bound")
  if(args[["hwe.threshold"]]>0.5 | args[["hwe.threshold"]] < 0) stop("hwe.threshold out off conventional bound")
  if(args[["colinearity.rsq"]]<0.8 | args[["colinearity.rsq"]]> 1) stop("colinearity.rsq out off conventional bound")
  if(!(args[["outputfile.format"]] %in% c("long","wide","specific"))) stop("invalid output.file.format");
  if(!(args[["output.exclude.covar"]] %in% c(0,1))) stop("output.exclude.covar must be 0 or 1");
  if(!(args[["Pvalue.type"]] %in% c("t.dist", "norm.dist"))) stop("Pvalue.type must be t.dist or norm.dist");
}




