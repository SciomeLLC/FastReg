#' fastR_hdf5convert a function to convert textual data to hdf5 format supported by FastReg()
#' @param data.file file name of the text data
#' @param h5.file file name of the result
#' @param header.row row containing the header
#' @param id.col column with the ids
#' @param data.col start of the data columns in the text data
#' @param buff.size buffer size
#' @param transpose boolean to transpose the data
#' @param chunk.edge chunk size for the hdf5 file
#' @return Boolean stating if the result file was created succesfully
#' @export
fastR_hdf5convert <- function(data.file, h5.file, header.row=1, id.col=1, data.col=2, buff.size=1.0, transpose=FALSE, chunk.edge=100) {
  if(!is.character(data.file) || length(data.file)>1) {
    cat("Error: data.file must be a single-member character vector\n")
    return(FALSE)
  }
  if(!is.character(h5.file) || length(h5.file)>1) {
    cat("Error: h5.file must be a single-member character vector\n")
    return(FALSE)
  }
  if(!is.numeric(header.row) || header.row!=as.integer(header.row) || length(header.row)>1 || header.row<1) {
    cat("Error: header.row must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  header.row <- as.integer(header.row)
  if(!is.numeric(id.col) || id.col!=as.integer(id.col) || length(id.col)>1 || id.col<1) {
    cat("Error: id.col must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  id.col <- as.integer(id.col)
  if(!is.numeric(data.col) || data.col!=as.integer(data.col) || length(data.col)>1 || data.col<=id.col) {
    cat("Error: header.row must be a single-member integer vector with value > id.col\n")
    return(FALSE)
  }
  data.col <- as.integer(data.col)
  if(!is.numeric(buff.size) || length(buff.size)>1 || buff.size<=0) {
    cat("Error: buff.size must be a single-member numeric vector with value > 0\n")
    return(FALSE)
  }
  buff.size <- as.double(buff.size)
  if(!is.logical(transpose) || length(transpose)>1) {
    cat("Error: transpose must be a single-member logical vector\n")
    return(FALSE)
  }
  if(!is.numeric(chunk.edge) || chunk.edge!=as.integer(chunk.edge) || length(chunk.edge)>1 || chunk.edge<1) {
    cat("Error: chunk.edge must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  chunk.edge <- as.integer(chunk.edge)
  return(.Call("hdf5convert", data.file, h5.file, header.row, id.col, data.col, buff.size, transpose, chunk.edge))
}
