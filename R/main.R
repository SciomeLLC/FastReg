#' FastReg a function to perform fast simple linear or logistic regression
#' @param config.file an optional character denoting configuration filename (default=NULL)
#' @param ... additional parameters not specified in config.file
#' @return numeric denoting elapsed.time
#' @export


FastReg <- function(config.file=NULL, ...) {
  protected.args <- list(...);
  # parse configuration file if exists
  if(!is.null(config.file)) args <- parseConfigFile(file=config.file, delim="\t", comment.char="#", protected.args=protected.args);
  # Assign default values if not specified in configuration file
  args <- assign.default.values(args=args);
  if(args[["verbose"]]>0) cat("Completed configuration parsing\n")

  # Set and combine covariate information
  names(args[["covariate.type"]]) <- names(args[["covariate.levels"]]) <- names(args[["covariate.ref.level"]]) <- names(args[["covariate.standardize"]]) <- args[["covariates"]]

  # validate configuration file and args
  validate.args(args);
  if(args[["verbose"]]>0) cat("Completed configuration preliminary validation\n")

  # load phenotype file as data frame
  pheno.df <-  as.data.frame(fread(file=args[["pheno.file"]], sep=args[["pheno.file.delim"]], header=TRUE),as.is=TRUE);
  if(!all(args[["pheno.rowname.cols"]]  %in% colnames(pheno.df))) stop("invalid pheno.row.name.cols");
  if(!(args[["phenotype"]] %in% colnames(pheno.df))) stop("invalid phenotype");
  pheno.df.rownames <- do.call(paste, c(lapply(args[["pheno.rowname.cols"]], function(x) pheno.df[[x]]), c("sep"="_")));

  # Get duplicates from phenotypes and subjects and filter them
  duplicate.subjects <- which(duplicated(pheno.df.rownames));
  num.duplicate.subjects <- length(duplicate.subjects);
  if(num.duplicate.subjects>0) {
    warning(paste(num.duplicate.subjects, "duplicate subject records were discared from pheno.file"));
    pheno.df.rownames <- pheno.df.rownames[-duplicate.subjects];
    pheno.df <- pheno.df[-duplicate.subjects,,drop=FALSE];
  }
  rownames(pheno.df) <- pheno.df.rownames;
  if(args[["verbose"]]>0) {
    cat("Completed pheno.file loading\n")
    cat("Discovered", nrow(pheno.df), "non-dulicate subjects in pheno.file\n");
  }
  remove(list=c("pheno.df.rownames", "num.duplicate.subjects", "duplicate.subjects"));

  # load covariates file as data frame
  covar.df <- as.data.frame(fread(file=args[["covar.file"]], sep=args[["covar.file.delim"]], header=TRUE),as.is=TRUE);
  if(!all(args[["covar.rowname.cols"]] %in% colnames(covar.df))) stop("invaild covar.rowname.cols");
  covar.df.rownames <- do.call(paste, c(lapply(args[["covar.rowname.cols"]], function(x) covar.df[[x]]), c("sep"="_")));
  # Get duplicates from covariates and subjects and filter them
  duplicate.subjects <- which(duplicated(covar.df.rownames));
  num.duplicate.subjects <- length(duplicate.subjects);
  if(num.duplicate.subjects>0) {
    warning(paste(num.duplicate.subjects, "duplicate subject records were discared from covar.file"));
    covar.df.rownames <- covar.df.rownames[-duplicate.subjects];
    covar.df <- covar.df[-duplicate.subjects,,drop=FALSE];
  }
  rownames(covar.df) <- covar.df.rownames;

  if(args[["verbose"]]>0) {
    cat("Completed covar.file loading\n")
    cat("Discovered", nrow(covar.df), "non-duplicate subjects in covar.file\n");
  }
  remove(list=c("covar.df.rownames", "num.duplicate.subjects", "duplicate.subjects"));



   ## convert txt POI file into h5 file if needed
  if(args[["POI.file.format"]]=="txt") {
    tmpfile <- gsub("[.]txt$", ".h5", args[["POI.file"]]);
    if(args[["verbose"]]>0) cat("Convering POI.file into h5 format\n")

    conversion.status <- try(convertPOI.TxtToH5(file=args[["POI.file.format"]], outfile=tmpfile, sep=args["POI.file.delim"]));
    if(inherits(conversion.status, "try-error")) stop("POI file TxtToH5 conversion failed")
    args[["POI.file"]] <- tmpfile;
    args[["POI.file.format"]] <- "h5";
  }

  # Get overlapping individuals and filter them
  POI.Individuals <- get.POI.individuals(type=args[["POI.file.format"]], file=args[["POI.file"]]);
  if(args[["verbose"]]>0) {
    cat("Discovered", length(POI.Individuals), "non-duplicate subjects in POI.file\n");
  }

  common.ind <- intersect(intersect(rownames(pheno.df), rownames(covar.df)), POI.Individuals);
  if(args[["verbose"]]>0) {
    cat(length(common.ind), "unique subjects where found to be common in pheno.file, covar.file and POI.file\n");
  }

  if(!is.null(args[["subject.subset.file"]])) {
    if(args[["verbose"]]>0) cat("Parsing subject.subset.file\n")

    subject.subset.df <-  as.data.frame(fread(file=args[["subject.subset.file"]], sep=args[["subject.subset.file.delim"]], header=TRUE),as.is=TRUE);

    if(!all(args[["subject.subset.rowname.cols"]]  %in% colnames(subject.subset.df))) stop("invalid subject.subset.row.name.cols");
    subset.Individuals <- do.call(paste, c(lapply(args[["subject.subset.rowname.cols"]], function(x) subject.subset.df[[x]]), c("sep"="_")));
    common.ind <- common.ind[common.ind %in% subset.Individuals];
    if(args[["verbose"]]>0) {
      cat(length(common.ind), "unique subjects where retained via subject.subset.file\n");
    }
    remove(list=c("subject.subset.df","subset.Individuals"))
  }


  if(length(common.ind)==0) stop("No overlapping individuals found in POI, pheno and covar files");
  POI.names <- get.POI.names(type=args[["POI.file.format"]], file=args[["POI.file"]]);
  if(args[["verbose"]]>0) {
    cat("Discovered", length(POI.names), "unique POIs in POI.file\n");
  }

  # Filter POI according to subset file if specified in configuration
  subset.POI <- POI.names;
  if(!is.null(args[["POI.subset.file"]])) {
    if(args[["verbose"]]>0) cat("Parsing POI.subset.file\n")

    POI.subset.df <-  as.data.frame(fread(file=args[["POI.subset.file"]], sep=args[["POI.subset.file.delim"]], header=TRUE),as.is=TRUE);

    if(!all(args[["POI.subset.rowname.cols"]]  %in% colnames(POI.subset.df))) stop("invalid POI.subset.row.name.cols");
    subset.POI <- do.call(paste, c(lapply(args[["POI.subset.rowname.cols"]], function(x) POI.subset.df[[x]]), c("sep"="_")));
    subset.POI <- POI.names[POI.names %in% subset.POI];
    if(args[["verbose"]]>0) {
      cat(length(subset.POI), "unique POIs where retained via POI.subset.file\n");
    }
    remove(list="POI.subset.df");
  }

  num.POI <- length(subset.POI);
  if(num.POI==0) stop("No overlapping individuals found in POI, pheno, covar files");



  POI.object <- open.POI.file(file=args[["POI.file"]], type=args[["POI.file.format"]], n=length(POI.Individuals), p=num.POI);
  poi.parser <- get.POI.Matrix(type=args[["POI.file.format"]], file.object=POI.object, POI.names=POI.names, POI.individuals=POI.Individuals);


  if(!is.null(args[["Split.by"]])) {
    if(!all(args[["Split.by"]] %in% colnames(covar.df))) stop("invalid Split.by argument");

    nsplit.vars <- length(args[["Split.by"]]);

    strata <- as.data.frame(do.call(table, covar.df[common.ind,args[["Split.by"]],drop=FALSE]), stringsAsFactors=FALSE, check.names=FALSE);
    colnames(strata) <- c(args[["Split.by"]], "Freq");
    strata <- strata[strata$Freq>0,];
    strata <- strata[do.call(order, strata[, args[["Split.by"]],drop=FALSE]),];
    nstrata <- nrow(strata);
    stratum.index.list <- list();
    stratum.id <- character(nstrata);
    for(stratum in 1:nstrata) {
        stratum.id[stratum] <- paste0("_",paste0(sapply(args[["Split.by"]], FUN=function(x) paste0(x,"=", strata[stratum, x])), collapse="_"));
        stratum.index.list[[stratum.id[stratum]]] <- common.ind[rowSums(do.call(cbind, lapply(args[["Split.by"]], FUN=function(x) covar.df[common.ind,x]==strata[stratum,x])))==nsplit.vars]
    }
  } else {
      nstrata <- 1;
      stratum.index.list <- list(common.ind);
      stratum.id <- "";
  }

 for(stratum in 1:nstrata) {

   outfile.Suffix <- stratum.id[stratum];

   if(!is.null(args[["Split.by"]]) & args[["verbose"]]>0) cat("Processing stratum:", sub("^_", "",stratum.id[stratum]), "\n");


  ind.set <- stratum.index.list[[stratum]];

  Y <- as.matrix(pheno.df[ind.set,args[["phenotype"]],drop=FALSE]);

  ## construct design matrix that includes covariate effects
  if(args[["verbose"]]>0) cat("Constructing Design Matrix\n")

  X <- createDesignMatrix(df=covar.df[ind.set,,drop=FALSE], covariates=args[["covariates"]],
                        covariate.type=args[["covariate.type"]],
                        covariate.standardize=args[["covariate.standardize"]],
                        no.intercept=args[["no.intercept"]],
                        covariate.levels=args[["covariate.levels"]],
                        covariate.ref.level=args[["covariate.ref.level"]],
                        colinearity.rsq=args[["colinearity.rsq"]], verbose=args[["verbose"]]);

  ### construct aux design matrix that includes main and interaction POI effects
  Z <- matrix(1, nrow=nrow(X), ncol=1, dimnames=list(ind.set, "Intercept"));

  if(!is.null(args[["POI.Covar.interactions"]])) {
    cZ <- do.call(c, lapply(args[["POI.Covar.interactions"]], function(x) which(regexpr(x,colnames(X))>0)));
    if(length(cZ)>0) Z <-cbind(Z, X[,cZ,drop=FALSE]);
    remove(list=c("cZ"));
  }

  # Get optimal block size if not specified in configuration. Checks for cores, multi-threading and memory and finds optimal block size
  if(args[["poi.block.size"]]==0) {
    if(args[["verbose"]]>0) cat("Estimating POI block size\n")
    poi.block.size <- estimate.poi.block.size(num.poi=num.POI, num.ind=length(ind.set), poi.type=args[["poi.type"]], num.cores=args[["num.cores"]])
  } else {
      poi.block.size <- args[["poi.block.size"]]
  }

  num.poi.blocks <- ceiling(num.POI/poi.block.size);
  if(args[["verbose"]]>0) {
    cat("POIs will be processed in ", num.poi.blocks, "blocks each of size", poi.block.size, "\n");
  }

  elapsed.time <- 0;
  # Output regression results to specified output file in chunks of specified block size
   for(block in 1:num.poi.blocks) {
     if(args[["verbose"]]>1) cat("Processing POIs block: ", block, "\n");

      poi.block <- subset.POI[(poi.block.size*(block-1)+1):min(poi.block.size*block, num.POI)];
      G <- poi.parser(poi=poi.block, ind=ind.set);

      if(args[["POI.type"]]=="genotypes") {
        if(args[["verbose"]]>1) cat("Performing MAF and HWE filtering\n");

        filter.res <- POI.filter(G=G, maf.threshold=args[["maf.threshold"]], hwe.threshold=args[["hwe.threshold"]]);

        fwrite(filter.res,  file=file.path(args[["output.dir"]], paste0("POI-Summary", outfile.Suffix,".tsv")),
                  sep="\t", append=(block>1));

        if(all(filter.res$keep==0)) {
          if(args[["verbose"]]>1) cat("no POI passed filtering\n");
          next;
        }
        G <- G[,filter.res$keep==1,drop=FALSE];
        filter.res <- filter.res[filter.res$keep==1, setdiff(colnames(filter.res), "keep")];


        if(args[["verbose"]]>1) cat("Performing POI transformation\n");
        G <- POI.transform(G, effect.type=args[["POI.effect.type"]]);
      } else {
        filter.res <- data.frame("POI"=poi.block,  stringsAsFactors=FALSE);
      }


      if(args[["regression.type"]]=="logistic") {
        if(args[["verbose"]]>1) cat("Performing logistic regression\n");
        res <- logisticRegression(Y=Y, G=G, X=X, Z=Z, max.iter=args[["max.iter"]]);
      } else {
        if(args[["verbose"]]>1) cat("Performing linear regression\n");
        res <- linearRegression(Y=Y, G=G, X=X, Z=Z);
      }

      nP <- nrow(res$beta);
      nR <- nrow(filter.res);

      res$pvl <- matrix(NA, nrow=nrow(res$beta), ncol=ncol(res$beta), dimnames=dimnames(res$beta));

      if(args[["Pvalue.type"]]=="t.dist") {
        res$pvl[] <- 2*(1-pt(abs(c(res$beta)/c(res$se_beta)), df=rep(res$DF, each=nP)));
      } else if(args[["Pvalue.type"]]=="norm.dist") {
        res$pvl[] <- 2*(1-pnorm(abs(c(res$beta)/c(res$se_beta))));
      }

      res$pvl[res$pvl>1] <- 1;
      elapsed.time <- elapsed.time + res$elapsed.time;


    if(args[["output.exclude.covar"]]==1) {
      rN <- which(regexpr("POI$",rownames(res[["beta"]]))>0);
      res[c("beta","se_beta", "pvl")] <- lapply(res[c("beta","se_beta")], function(x) x[rN,,drop=FALSE]);
    }


    if(args[["verbose"]]>1) cat("saving results\n");

    if(args[["output.file.format"]]=="long") {
      res <- cbind(filter.res[rep(1:nR, each=nP), "POI",drop=FALSE],  N=rep(res$N, each=nP), DF=rep(res$DF, each=nP), "Effect"=rep(rownames(res$beta), nR), "Estimate"=c(res$beta), "Std Error"=c(res$se_beta), "P-value"=c(res$pvl));
      fwrite(res, file=file.path(args[["output.dir"]], paste0("Results", outfile.Suffix,".tsv")), append=(block>1), sep="\t");
    } else if(args[["output.file.format"]]=="wide") {
      rownames(res$beta) <- paste(rownames(res$beta), "Estimate");
      rownames(res$se_beta) <- paste(rownames(res$se_beta), "Std Error");
      rownames(res$pvl) <- paste(rownames(res$pvl), "P-value");

      res <- cbind(filter.res[, "POI",drop=FALSE],  N=res$N, DF=res$DF, "Estimate"=t(res$beta), "Std Error"=t(res$se_beta), "P-value"=t(res$pvl));
      fwrite(res, file=file.path(args[["output.dir"]], paste0("Results", outfile.Suffix, ".tsv")), append=(block>1), sep="\t");
    } else {
        for(en in colnames(res$beta)) {
          fwrite(cbind(filter.res[, "POI", drop=FALSE],  N=res$N, DF=res$DF, "Estimate"=res$beta[en,], "Std Error"=res$se_beta[en,], "P-value"=res$pvl[en,]),
                    file=file.path(args[["outdir"]], paste0(en, outfile.Suffix, ".tsv")), sep="\t", append=(block>1));
        }
    }

   }
  }

  ## close POI file connection
  close.POI.file(type=args[["POI.file.format"]], POI.object);

  invisible(elapsed.time)

}
