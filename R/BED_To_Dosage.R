#' Bed_To_dosage function to convert previously generated H5 files into plink files
#' @param input.dir folder location where input files reside
#' @param output.dir folder location where plink files will be saved
#' @param prefix a character denoting which set of h5 files to convert (first set is used by default)
#' @param poi.chunk.size integer denoting chunk.size to be used during conversion
#' @return integer 1 when successful
#' @import BEDMatrix
#' @import data.table
library("BEDMatrix");
library("data.table");
Bed_To_dosage <- function(input.dir = ".", output.dir = NULL, prefix = NULL, poi.chunk.size = 1000) {
  
  if(is.null(output.dir)) output.dir <- input.dir;
  
  available.files <- list.files(path=input.dir, pattern="[.]bed");
  
  if(length(available.files)==0) stop("input.dir does not consist of any .bed file");
  if(is.null(prefix)) prefix <- sub("[.]bed", "", basename(available.files)[1]);
  
  if(length(prefix)>1) {
      warning("only first prefix used")
      prefix <- prefix[1];
  }
  
  src.files <- file.path(input.dir, paste0(prefix, ".", c("bim", "fam", "bed")));
  out.files <- file.path(output.dir, paste0(prefix, ".", c("dosage.txt")));
  
  num.poi <- as.integer(strsplit(system(command=paste("wc -l", src.files[1]), intern=TRUE),split=" ")[[1]][1]);
  num.ind <- as.integer(strsplit(system(command=paste("wc -l", src.files[2]), intern=TRUE),split=" ")[[1]][1]);
  
  num.poi.chunks <- ceiling(num.poi/poi.chunk.size);
  poi.chunk.size <- ceiling(num.poi/num.poi.chunks);
  
  
  bm <- BEDMatrix(src.files[3], n=num.ind, p=num.poi);
  
  header <- c("SNP", "A1","A2", paste0(paste0("FID",formatC(1:num.ind, format="d", flag="0", digits=floor(log10(num.ind)))), "\t",paste0("IND", formatC(1:num.ind, format="d", flag="0", digits=floor(log10(num.ind))))));
  
  
  for( i in 1:num.poi.chunks) {
		start.pos <- (i-1)*poi.chunk.size+1;
		end.pos <- min(i*poi.chunk.size, num.poi);
		nr <- end.pos-start.pos+1;
		df <- do.call(rbind, strsplit(system(command=paste0("sed -n '", start.pos,",",end.pos,"p' ", src.files[1]), intern=TRUE),split="\t"));
		df <- as.data.frame(df, stringsAsFactors=FALSE, check.names=FALSE);
		# unlink(tmpfile);
		colnames(df) <- c("seq", "id", "strand", "pos", "ref", "alt");
		poi.block <- t(bm[,seq(start.pos, end.pos),drop=FALSE]);
		val <- c("0"="1.00\t0.00", "1"="0.00\t1.00","2"="0.00\t0.00");
		dosage.block <- matrix("NA\tNA", nrow=nr, ncol=num.ind);
		dosage.block[!is.na(poi.block)] <- val[as.character(poi.block[!is.na(poi.block)])];
		dosage.block <- cbind(df[,c("id", "ref", "alt")], dosage.block);
    print(dim(dosage.block))
    print(length(header))
		colnames(dosage.block) <- header;
		fwrite(dosage.block, file=out.files[1], row.names=FALSE, col.names=(i==1), sep="\t", append=(i>1), quote=FALSE); 	
  }
  
  # close(con);
  
  invisible(1);
}
