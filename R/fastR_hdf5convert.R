#' hdf5convert a function to convert textual data to hdf5 format supported by FastReg()
#' @param dataFile file name of the text data
#' @param h5Dir file dir to write hdf5 files to
#' @param headerRow row containing the header
#' @param idCol column with the ids
#' @param dataCol start of the data columns in the text data
#' @param buffSize buffer size
#' @param transpose boolean to transpose the data
#' @param chunkEdge chunk size for the hdf5 file
#' @param vcf TRUE or FALSE
#' @param delimiter delimiter used in the dataFile
#' @param gz bool to gzip the h5 dataset or not
#' @param poiPerFile =-1
#' @param singleFile =FALSE
#' @param serverThreads =-1
#' @param serverMem =-1.0
#' @return Boolean stating if the result file was created succesfully
#' @export

fastR_hdf5convert <- function(
    dataFile,
    h5Dir,
    headerRow = 1,
    idCol = 1,
    dataCol = 2,
    buffSize = 1.0,
    transpose = FALSE,
    chunkEdge = 100,
    vcf = FALSE,
    delimiter = " \t",
    gz = FALSE,
    poiPerFile = -1,
    singleFile = FALSE,
    serverThreads = -1,
    serverMem = -1.0) {
  if (!is.character(dataFile) || length(dataFile) > 1) {
    cat("Error: dataFile must be a single-member character vector\n")
    return(FALSE)
  }
  if (!is.character(h5Dir) || length(h5Dir) > 1) {
    cat("Error: h5Dir must be a single-member character vector\n")
    return(FALSE)
  } else if (file.exists(h5Dir)) {
    cat("Error: file or folder named \"", h5Dir, "\" already exists\n", sep = "")
    return(FALSE)
  }
  dir.create(h5Dir)
  if (!is.numeric(headerRow) || headerRow != as.integer(headerRow) || length(headerRow) > 1 || headerRow < 1) {
    cat("Error: headerRow must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  headerRow <- as.integer(headerRow)
  if (!is.numeric(idCol) || idCol != as.integer(idCol) || length(idCol) > 1 || idCol < 1) {
    cat("Error: idCol must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  idCol <- as.integer(idCol)
  if (!is.numeric(dataCol) || dataCol != as.integer(dataCol) || length(dataCol) > 1 || dataCol <= idCol) {
    cat("Error: headerRow must be a single-member integer vector with value > idCol\n")
    return(FALSE)
  }
  dataCol <- as.integer(dataCol)
  if (!is.numeric(buffSize) || length(buffSize) > 1 || buffSize <= 0) {
    cat("Error: buffSize must be a single-member numeric vector with value > 0\n")
    return(FALSE)
  }
  buffSize <- as.double(buffSize)
  if (!is.logical(transpose) || length(transpose) > 1) {
    cat("Error: transpose must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(chunkEdge) || chunkEdge != as.integer(chunkEdge) || length(chunkEdge) > 1 || chunkEdge < 1) {
    cat("Error: chunkEdge must be a single-member integer vector with value > 0\n")
    return(FALSE)
  }
  chunkEdge <- as.integer(chunkEdge)
  if (!is.logical(vcf) || length(vcf) > 1) {
    cat("Error: vcf must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.character(delimiter) || length(delimiter) > 1) {
    cat("Error: delimiter must be a single-member character vector\n")
    return(FALSE)
  }
  if (!is.logical(gz) || length(gz) > 1) {
    cat("Error: gz must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(poiPerFile) || poiPerFile != as.integer(poiPerFile) || length(poiPerFile) > 1 || poiPerFile < chunkEdge) {
    if (poiPerFile != -1) {
      cat("Error: poiPerFile must be a single-member integer vector with value >= chunkEdge\n")
      return(FALSE)
    }
  }
  poiPerFile <- as.integer(poiPerFile)
  if (!is.logical(singleFile) || length(singleFile) > 1) {
    cat("Error: singleFile must be a single-member logical vector\n")
    return(FALSE)
  }
  if (!is.numeric(serverThreads) || serverThreads != as.integer(serverThreads) || length(serverThreads) > 1) {
    cat("Error: serverThreads must be a single-member integer vector\n")
    return(FALSE)
  }
  serverThreads <- as.integer(serverThreads)
  if (!is.numeric(serverMem) || length(serverMem) > 1) {
    cat("Error: serverMem must be a single-member numeric vector\n")
    return(FALSE)
  }
  serverMem <- as.double(serverMem)
  return(.Call(
    "fastR_hdf5convert_cpp",
    h5Dir,
    headerRow,
    idCol,
    dataCol,
    buffSize,
    transpose,
    chunkEdge,
    vcf,
    delimiter,
    gz,
    poiPerFile,
    singleFile,
    serverThreads,
    serverMem
  ))
}
