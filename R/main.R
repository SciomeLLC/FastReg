#' FastReg a function to perform fast simple linear or logistic regression
#' @param config.file an optional character denoting configuration filename (default=NULL)
#' @param ... additional parameters not specified in config.file
#' @return numeric denoting elapsed.time
#' @export

FastReg <- function(config.file=NULL, ...) {
  protected.args <- list(...);

  if(!is.null(config.file)) args <- parseConfigFile(file=config.file, delim="\t", comment.char="#", protected.args=protected.args);

  args <- assign.default.values(args=args);
  if(args[["verbose"]]>0) cat("Completed configuration parsing\n")


  names(args[["covariate.type"]]) <- names(args[["covariate.levels"]]) <- names(args[["covariate.ref.level"]]) <- names(args[["covariate.standardize"]]) <- args[["covariates"]]

  validate.args(args);
  if(args[["verbose"]]>0) cat("Completed configuration preliminary validation\n")

  pheno.df <-  as.data.frame(fread(file=args[["pheno.file"]], sep=args[["pheno.file.delim"]], header=TRUE),as.is=TRUE);
  if(!all(args[["pheno.rowname.cols"]]  %in% colnames(pheno.df))) stop("invalid pheno.row.name.cols");
  if(!(args[["phenotype"]] %in% colnames(pheno.df))) stop("invalid phenotype");
  pheno.df.rownames <- do.call(paste, c(lapply(args[["pheno.rowname.cols"]], function(x) pheno.df[[x]]), c("sep"="_")));
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

  covar.df <- as.data.frame(fread(file=args[["covar.file"]], sep=args[["covar.file.delim"]], header=TRUE),as.is=TRUE);
  if(!all(args[["covar.rowname.cols"]] %in% colnames(covar.df))) stop("invaild covar.rowname.cols");
  covar.df.rownames <- do.call(paste, c(lapply(args[["covar.rowname.cols"]], function(x) covar.df[[x]]), c("sep"="_")));
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

  pheno.df <- pheno.df[common.ind,,drop=FALSE];
  covar.df <- covar.df[common.ind,,drop=FALSE];
  Y <- as.matrix(pheno.df[,args[["phenotype"]],drop=FALSE]);

  ## construct design matrix that includes covariate effects
  if(args[["verbose"]]>0) cat("Constructing Design Matrix\n")

  X <- createDesignMatrix(df=covar.df, covariates=args[["covariates"]],
                        covariate.type=args[["covariate.type"]],
                        covariate.standardize=args[["covariate.standardize"]],
                        no.intercept=args[["no.intercept"]],
                        covariate.levels=args[["covariate.levels"]],
                        covariate.ref.level=args[["covariate.ref.level"]],
                        colinearity.rsq=args[["colinearity.rsq"]], verbose=args[["verbose"]]);

  ### construct aux design matrix that includes main and interaction POI effects
  Z <- matrix(1, nrow=nrow(X), ncol=1, dimnames=list(common.ind, "Intercept"));

  if(!is.null(args[["POI.Covar.interactions"]])) {
    cZ <- do.call(c, lapply(args[["POI.Covar.interactions"]], function(x) which(regexpr(x,colnames(X))>0)));
    if(length(cZ)>0) Z <-cbind(Z, X[,cZ,drop=FALSE]);
    remove(list=c("cZ"));
  }

  if(args[["poi.block.size"]]==0) {
    ### place holder function to compute poi.block.size
    if(args[["verbose"]]>0) cat("Estimating POI block size\n")

    args[["poi.block.size"]] <- estimate.poi.block.size(num.poi=num.POI, num.ind=length(common.ind), poi.type=args[["poi.type"]], num.cores=args[["num.cores"]])
  }

  num.poi.blocks <- ceiling(num.POI/args[["poi.block.size"]]);
  if(args[["verbose"]]>0) {
    cat("POIs will be processed in ", num.poi.blocks, "blocks each of size", args[["poi.block.size"]], "\n");
  }
  poi.parser <- get.POI.Matrix(type=args[["POI.file.format"]], file.object=POI.object, POI.names=POI.names, POI.individuals=POI.Individuals);

  elapsed.time <- 0;

   for(block in 1:num.poi.blocks) {
     if(args[["verbose"]]>1) cat("Processing POIs block: ", block, "\n");

      poi.block <- subset.POI[(args[["poi.block.size"]]*(block-1)+1):min(args[["poi.block.size"]]*block, num.POI)];
      G <- poi.parser(poi=poi.block, ind=common.ind);

      if(args[["POI.type"]]=="genotypes") {
        if(args[["verbose"]]>1) cat("Performing MAF and HWE filtering\n");

        filter.res <- POI.filter(G=G, maf.threshold=args[["maf.threshold"]], hwe.threshold=args[["hwe.threshold"]]);

        fwrite(filter.res[filter.res$keep==0, setdiff(colnames(filter.res), "keep")], file=file.path(args[["output.dir"]], "Filtered-POI.tsv"), sep="\t", append=(block>1));
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

      elapsed.time <- elapsed.time + res$elapsed.time;

    nP <- ncol(res$beta);
    nR <- nrow(filter.res);

    if(args[["output.exclude.covar"]]==1) {
      rN <- which(regexpr("POI$",colnames(res[["beta"]]))>0);
      res[c("beta","se_beta", "pvl")] <- lapply(res[c("beta","se_beta","pvl")], function(x) x[,rN,drop=FALSE]);
    }


    if(args[["verbose"]]>1) cat("saving results\n");

    if(args[["output.file.format"]]=="long") {
      res <- cbind(filter.res[rep(1:nR, each=ncol(res$beta)),],  N=rep(res$N, each=nP), DF=rep(res$DF, each=nP), "Effect"=rep(colnames(res$beta), nR), "Estimate"=c(t(res$beta)), "Std Error"=c(t(res$se_beta)), "P-value"=c(t(res$pvl)));
      fwrite(res, file=file.path(args[["output.dir"]], "Results.tsv"), append=(block>1), sep="\t");
    } else if(args[["output.file.format"]]=="wide") {
      colnames(res$beta) <- paste(colnames(res$beta), "Estimate");
      colnames(res$se_beta) <- paste(colnames(res$se_beta), "Std Error");
      colnames(res$pvl) <- paste(colnames(res$pvl), "P-value");

      res <- cbind(filter.res,  N=res$N, DF=res$DF, "Estimate"=res$beta, "Std Error"=res$se_beta, "P-value"=res$pvl);
      fwrite(res, file=file.path(args[["output.dir"]], "Results.tsv"), append=(block>1), sep="\t");
    } else {
        for(en in colnames(res$beta)) {
          fwrite(cbind(filter.res,  N=res$N, DF=res$DF, "Estimate"=res$beta[,en], "Std Error"=res$se_beta[,en], "P-value"=res$pvl[,en]),
                    file=file.path(args[["outdir"]], paste0(en, ".tsv")), sep="\t", append=(block>1));
        }
    }

   }

  ## close POI file connection
  close.POI.file(type=args[["POI.file.format"]], POI.object);

  invisible(elapsed.time)

}
