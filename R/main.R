#' FastReg a function to perform fast simple linear or logistic regression
#' @param config.file an optional character denoting configuration filename (default=NULL)
#' @param ... additional parameters not specified in config.file
#' @return numeric denoting elapsed.time
#' @export

# 8 digit precision for results
FastReg <- function(config.file=NULL, ...) {
  protected.args <- list(...);

  if(!is.null(config.file)) args <- parseConfigFile(file=config.file, delim="\t", comment.char="#", protected.args=protected.args);

  args <- assign.default.values(args=args);
  names(args[["covariate.type"]]) <- names(args[["covariate.levels"]]) <- names(args[["covariate.ref.level"]]) <- names(args[["covariate.standardize"]]) <- args[["covariates"]]

  validate.args(args);

  pheno.df <-  as.data.frame(fread(file=args[["pheno.file"]], sep=args[["pheno.file.delim"]], header=TRUE),as.is=TRUE);
  if(!all(args[["pheno.rowname.cols"]]  %in% colnames(pheno.df))) stop("invalid pheno.row.name.cols");
  if(!(args[["phenotype"]] %in% colnames(pheno.df))) stop("invalid phenotype");
  rownames(pheno.df) <- do.call(paste, c(lapply(args[["pheno.rowname.cols"]], function(x) pheno.df[[x]]), c("sep"="_")));


  covar.df <- as.data.frame(fread(file=args[["covar.file"]], sep=args[["covar.file.delim"]], header=TRUE),as.is=TRUE);
  if(!all(args[["covar.rowname.cols"]] %in% colnames(covar.df))) stop("invaild covar.rowname.cols");
  rownames(covar.df) <- do.call(paste, c(lapply(args[["covar.rowname.cols"]], function(x) covar.df[[x]]), c("sep"="_")));


   ## convert txt POI file into h5 file if needed
  if(args[["POI.file.format"]]=="txt") {
    tmpfile <- gsub("[.]txt$", ".h5", args[["POI.file"]]);
    conversion.status <- try(convertPOI.TxtToH5(file=args[["POI.file.format"]], outfile=tmpfile, sep=args["POI.file.delim"]));
    if(inherits(conversion.status, "try-error")) stop("POI file TxtToH5 conversion failed")
    args[["POI.file"]] <- tmpfile;
    args[["POI.file.format"]] <- "h5";
  }

  POI.Individuals <- get.POI.individuals(type=args[["POI.file.format"]], file=args[["POI.file"]]);

  common.ind <- intersect(intersect(rownames(pheno.df), rownames(covar.df)), POI.Individuals);
  if(length(common.ind)==0) stop("No overlapping individuals found in POI, pheno, covar files");
  POI.names <- get.POI.names(type=args[["POI.file.format"]], file=args[["POI.file"]]);
  num.POI <- length(POI.names);

  POI.object <- open.POI.file(file=args[["POI.file"]], type=args[["POI.file.format"]], n=length(POI.Individuals), p=num.POI);

  pheno.df <- pheno.df[common.ind,,drop=FALSE];
  covar.df <- covar.df[common.ind,,drop=FALSE];
  Y <- as.matrix(pheno.df[,args[["phenotype"]],drop=FALSE]);
  ## construct design matrix that includes covariate effects
  X <- createDesignMatrix(df=covar.df, covariates=args[["covariates"]],
                        covariate.type=args[["covariate.type"]],
                        covariate.standardize=args[["covariate.standardize"]],
                        no.intercept=args[["no.intercept"]],
                        covariate.levels=args[["covariate.levels"]],
                        covariate.ref.level=args[["covariate.ref.level"]],
                        colinearity.rsq=args[["colinearity.rsq"]]);

  ### construct aux design matrix that includes main and interaction POI effects
  Z <- matrix(1, nrow=nrow(X), ncol=1, dimnames=list(common.ind, "Intercept"));

  if(!is.null(args[["POI.covar.interactions"]])) {
    cZ <- do.call(c, lapply(args[["POI.covar.interactions"]], which(regexpr(x,colnames(X))>0)));
    if(length(cZ)>0) Z <-cbind(Z, X[,cZ,drop=FALSE]);
    remove(list=c("cZ"));
  }

  if(args[["poi.block.size"]]==0) {
    ### place holder function to compute poi.block.size
    args[["poi.block.size"]] <- 50
    # args[["poi.block.size"]] <- estimate.poi.block.size(num.poi=num.POI, num.ind=length(common.ind), poi.type=args[["poi.type"]], num.cores=args[["num.cores"]])
  }

  num.poi.blocks <- ceiling(length(POI.names)/args[["poi.block.size"]]);

  poi.parser <- get.POI.Matrix(type=args[["POI.file.format"]], file.object=POI.object, POI.names=POI.names, POI.individuals=POI.Individuals);

  elapsed.time <- 0;
  poi.block.size <- args[["poi.block.size"]];
   for(block in 1:num.poi.blocks) {
      poi.block <- POI.names[(poi.block.size*(block-1)+1):min(poi.block.size*block, num.POI)];
      G <- poi.parser(poi=poi.block, ind=common.ind);
      filter.res <- POI.filter(G=G, maf.threshold=args[["maf.threshold"]], hwe.threshold=args[["hwe.threshold"]]);

      fwrite(filter.res[filter.res$keep==0,], file=file.path(args[["output.dir"]], "Filtered-POI.tsv"), sep="\t", append=(block>1));


      if(all(filter.res$keep==0)) next;
      G <- G[,filter.res$keep==1,drop=FALSE];

      G <- POI.transform(G, effect.type=args[["POI.effect.type"]]);

      if(args[["regression.type"]]=="logistic") {
        res <- logisticRegression(Y=Y, G=G, X=X, Z=Z);
      } else {
        res <- linearRegression(Y=Y, G=G, X=X, Z=Z);
      }

      elapsed.time <- elapsed.time + res$elapsed.time;

    nP <- ncol(res$beta);
    if(args[["output.exclude.covar"]]==1) {
      rN <- which(regexpr("POI$",colnames(res[["beta"]]))>0);
      res[c("beta","se_beta", "pvl")] <- lapply(res[c("beta","se_beta","pvl")], function(x) x[,rN,drop=FALSE]);
    }



    if(args[["output.file.format"]]=="long") {
      res <- cbind(filter.res[rep(which(filter.res$keep==1), each=ncol(res$beta)),],  N=rep(res$N, each=ncol(res$beta)), DF=rep(res$DF, each=ncol(res$beta)), "Effect"=rep(colnames(res$beta), nrow(res$beta)), "Estimate"=c(t(res$beta)), "Std Error"=c(t(res$se_beta)), "P-value"=c(t(res$pvl)));
      fwrite(res, file=file.path(args[["output.dir"]], "Results.tsv"), append=(block>1), sep="\t");
    } else if(args[["output.file.format"]]=="wide") {
      colnames(res$beta) <- paste(colnames(res$beta), "Estimate");
      colnames(res$se_beta) <- paste(colnames(res$se_beta), "Std Error");
      colnames(res$pvl) <- paste(colnames(res$pvl), "P-value");

      res <- cbind(filter.res[filter.res$keep==1,],  N=res$N, DF=res$DF, "Estimate"=res$beta, "Std Error"=res$se_beta, "P-value"=res$pvl);
      fwrite(res, file=file.path(args[["output.dir"]], "Results.tsv"), append=(block>1), sep="\t");
    } else {
        for(en in colnames(res$beta)) {
          fwrite(cbind(filter.res[filter.res$keep==1,],  N=res$N, DF=res$DF, "Estimate"=res$beta[,en], "Std Error"=res$se_beta[,en], "P-value"=res$pvl[,en]),
                    file=file.path(args[["outdir"]], paste0(en, ".tsv")), sep="\t", append=(block>1));
        }
    }

   }

  ## close POI file connection
  close.POI.file(type=args[["POI.file.format"]], POI.object);

  invisible(elapsed.time)

}
