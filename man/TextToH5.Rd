% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{TextToH5}
\alias{TextToH5}
\title{TextToH5 a function to convert textual data to hdf5 format supported by FastReg()}
\usage{
TextToH5(
  data.file,
  h5.dir,
  header.row = 1,
  id.col = 1,
  data.col = 2,
  buff.size = 1,
  transpose = FALSE,
  chunk.edge = 100,
  vcf = FALSE,
  delimiter = " \\t",
  gz = FALSE,
  poi.per.file = -1,
  single.file = FALSE,
  server.threads = -1,
  server.mem = -1
)
}
\arguments{
\item{data.file}{Path of the input text file.}

\item{h5.dir}{Path of a new directory to write generated hdf5 files.}

\item{header.row}{Default: 1. Row number of header line.}

\item{id.col}{Default 1. Column number that contains identifiers.}

\item{data.col}{Default: 2. Column number that contains first data entry.}

\item{buff.size}{Default: 1.0. Size in Gb of memory available for data conversion.}

\item{transpose}{Default: FALSE. Boolean value that indicates whether text data should be transposed. Set TRUE when columns represent individuals and rows represent POIs.}

\item{chunk.edge}{Default: 100. Specify size of data chunks to be applied to output hdf5 files. Must be less than or equal to values of poi.per.file and greater than zero.}

\item{vcf}{Default: FALSE. Indicate whether input file is in vcf format. If TRUE, header.row, id.col, data.col and transpose values are ignored.}

\item{delimiter}{Specify a string of all delimiters used in the text file to separate columns. Default will split entries by space or tab.}

\item{gz}{Default: FALSE. Indicate whether input text file is gzip-compressed.}

\item{poi.per.file}{Default: -1. Indicate the number of POIs to write to each output hdf5 file. A value of -1 indicates the count should be calculated based on available system resources.}

\item{single.file}{Default: FALSE. Indicate whether a single hdf5 file should be produced rather than a series.}

\item{server.threads}{Default: -1. Indicate the total number of CPU threads available for FastReg - utilized to determine the optimal number of output files. A value of -1 indicates the count should be detected automatically. Set this parameter when running FastReg on a shared system and specify resources allocated for your own use only.}

\item{server.mem}{Default: -1.0. Indicate the total memory in Gb available for FastReg - utilized to determine the optimal number of output files. A value of -1.0 indicates the total should be detected automatically. Set this parameter when running FastReg on a shared system and specify resources allocated for your own use only.}
}
\value{
0 = success, 1 = failure
}
\description{
TextToH5 a function to convert textual data to hdf5 format supported by FastReg()
}
