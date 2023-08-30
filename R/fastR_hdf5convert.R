#' hdf5convert a function to convert textual data to hdf5 format supported by FastReg()
#' @param dataFile file name of the text data
#' @param h5File file name of the result
#' @param headerRow row containing the header
#' @param idCol column with the ids
#' @param dataCol start of the data columns in the text data
#' @param buffSize buffer size
#' @param transpose boolean to transpose the data
#' @param chunkEdge chunk size for the hdf5 file
#' @return Boolean stating if the result file was created succesfully
#' @export
fastR_hdf5convert <- function(dataFile, h5File, headerRow=1, idCol=1, dataCol=2, buffSize=1.0, transpose=FALSE, chunkEdge=100) {
  if(!is.character(dataFile) || length(dataFile)>1) {
    cat("Error: dataFile must be a single-member character vector\n")
    return(FALSE)
  }
  if(!is.character(h5File) || length(h5File)>1) {
    cat("Error: h5File must be a single-member character vector\n")
    return(FALSE)
  }
  if(!is.numeric(headerRow) || headerRow!=as.integer(headerRow) || length(headerRow)>1 || headerRow<1) {
    cat("Error: headerRow must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  headerRow <- as.integer(headerRow)
  if(!is.numeric(idCol) || idCol!=as.integer(idCol) || length(idCol)>1 || idCol<1) {
    cat("Error: idCol must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  idCol <- as.integer(idCol)
  if(!is.numeric(dataCol) || dataCol!=as.integer(dataCol) || length(dataCol)>1 || dataCol<=idCol) {
    cat("Error: headerRow must be a single-member integer vector with value > idCol\n")
    return(FALSE)
  }
  dataCol <- as.integer(dataCol)
  if(!is.numeric(buffSize) || length(buffSize)>1 || buffSize<=0) {
    cat("Error: buffSize must be a single-member numeric vector with value > 0\n")
    return(FALSE)
  }
  buffSize <- as.double(buffSize)
  if(!is.logical(transpose) || length(transpose)>1) {
    cat("Error: transpose must be a single-member logical vector\n")
    return(FALSE)
  }
  if(!is.numeric(chunkEdge) || chunkEdge!=as.integer(chunkEdge) || length(chunkEdge)>1 || chunkEdge<1) {
    cat("Error: chunkEdge must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  chunkEdge <- as.integer(chunkEdge)
  return(.Call("hdf5convert", dataFile, h5File, headerRow, idCol, dataCol, buffSize, transpose, chunkEdge))
}
