#include <chunker.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <hdf5.h>
#include <zlib.h>

using namespace std;
using namespace Rcpp;

// consolidates information derived from data file preprocessing
struct preprocess_info
{
  char *collabel;
  char *rowlabel;
  char **colnames;
  char *colline;
  hsize_t dimensions[2];
  hsize_t row_dim;
  int growdim;
  int fixdim;
};

// consolidates all HDF5 handles and attributes
struct hdf5_vars
{
  hid_t file;
  hid_t vals_dataset;
  hid_t vals_datatype;
  hid_t vals_dataspace;
  hid_t vals_memspace;
  hid_t vals_prop_list;
  hid_t col_dataset;
  hid_t col_datatype;
  hid_t col_dataspace;
  hid_t col_prop_list;
  hid_t row_dataset;
  hid_t row_datatype;
  hid_t row_dataspace;
  hid_t row_memspace;
  hid_t row_prop_list;
  int written;
};

// consolidates per-HDF5 file dimension variables
struct dim_vars
{
  hsize_t vals_dataspace_dims[2];
  hsize_t vals_memspace_dims[2];
  hsize_t vals_hyperslab_pos[2];
  hsize_t vals_chunk_dims[2];
};

// consolidates user-defined parameters
struct fastR_user_params
{
  char *infile_name;
  char *h5file_base;
  size_t header_row;
  size_t name_column;
  size_t data_column;
  size_t data_buffer_max;
  int transpose;
  hsize_t chunk_edge;
  int vcf;
  char *delim;
  int gz;
  int poi_per_file;
  int single;
  int server_threads;
  float server_memory;
};

// consolidates read buffer pointers
struct read_buffers
{
  float ***val_buffer;
  char **row_buffer;
};

// implements POSIX function strsep for portability
char *separate(char **stringp, const char *delim)
{
  size_t span;
  char *first, *last;
  if (*stringp == NULL)
  {
    return (NULL);
  }
  first = *stringp;
  span = strcspn(first, delim);
  last = first + span;
  if (*last == '\0')
  {
    *stringp = NULL;
  }
  else
  {
    *last = '\0';
    *stringp = last + 1;
  }
  return (first);
}

// implements POSIX function getline for portability - adds gz read capability,
// skips empty lines
ssize_t get_full_line(char **lineptr, size_t *n, FILE *stream, int gzflag,
                      gzFile gzstream)
{
  char *res, *start, *mark;
  size_t prev = 0, count;
  if (*lineptr == NULL || *n == 0)
  {
    *lineptr = (char *)malloc(10485760);
    *n = 10485760;
    if (*lineptr == NULL)
    {
      return (-1);
    }
  }
  start = *lineptr;
  count = *n;
  if (gzflag == 1)
  {
    res = gzgets(gzstream, start, (int)count);
    while (res != NULL && start[0] == '\n')
    {
      res = gzgets(gzstream, start, (int)count);
    }
  }
  else
  {
    res = fgets(start, (int)count, stream);
    while (res != NULL && start[0] == '\n')
    {
      res = fgets(start, (int)count, stream);
    }
  }
  if (res == NULL)
  {
    return (-1);
  }
  else
  {
    while (1)
    {
      mark = (char *)memchr((void *)start, '\n', count);
      if (mark != NULL)
      {
        return ((ssize_t)(mark - *lineptr) + 1);
      }
      prev = *n;
      *n = 2 * prev;
      *lineptr = (char *)realloc(*lineptr, *n);
      if (*lineptr == NULL)
      {
        return (-1);
      }
      start = *lineptr + prev - 1;
      count = prev + 1;
      if (gzflag == 1)
      {
        res = gzgets(gzstream, start, (int)count);
      }
      else
      {
        res = fgets(start, (int)count, stream);
      }
      if (res == NULL)
      {
        return (-1);
      }
    }
  }
}

static int preprocess_datafile(FILE *datafile, gzFile gzdatafile,
                               struct fastR_user_params *par,
                               struct preprocess_info *pi)
{
  ssize_t nread;
  char *line = NULL;
  char *buff, *pch = NULL;
  size_t colnamessize;
  size_t len = 0;
  hsize_t i, j, k;
  i = 0;
  j = 0;
  k = 0;
  Rcpp::Rcout << "Pre-processing input file..." << std::endl;
  // set column data set name based on matrix orientation
  if (par->transpose == 1)
  {
    pi->collabel = strdup("individuals");
    pi->rowlabel = strdup("predictors_of_interest");
  }
  else
  {
    pi->rowlabel = strdup("individuals");
    pi->collabel = strdup("predictors_of_interest");
  }
  pi->colnames = (char **)malloc(sizeof(char *));
  colnamessize = 1;
  // toss lines until header row is read
  if (par->vcf == 0)
  {
    nread = get_full_line(&pi->colline, &len, datafile, par->gz, gzdatafile);
    i++;
    if (nread == -1)
    {
      Rprintf("Error: could not read specified header row from input file\n");
      free(line);
      return (1);
    }
    while (i < par->header_row)
    {
      free(pi->colline);
      pi->colline = NULL;
      len = 0;
      nread = get_full_line(&pi->colline, &len, datafile, par->gz, gzdatafile);
      i++;
      if (nread == -1)
      {
        Rprintf("Error: could not read specified header row from input file\n");
        free(line);
        return (1);
      }
    }
  }
  else
  {
    nread = get_full_line(&pi->colline, &len, datafile, par->gz, gzdatafile);
    i++;
    if (nread == -1)
    {
      Rprintf("Error: could not read specified header row from input file\n");
      free(line);
      return (1);
    }
    while (pi->colline[0] == '#' && pi->colline[1] == '#')
    {
      free(pi->colline);
      pi->colline = NULL;
      len = 0;
      nread = get_full_line(&pi->colline, &len, datafile, par->gz, gzdatafile);
      i++;
      if (nread == -1)
      {
        Rprintf("Error: could not read specified header row from input file\n");
        free(line);
        return (1);
      }
    }
    par->header_row = i;
  }
  // parse and store column names
  if (pi->colline[nread - 1] == '\n')
  {
    pi->colline[nread - 1] = '\0';
  }
  // hold line pointer so memory can be freed after parsing
  buff = pi->colline;
  // toss fields before first data point
  while (k < par->data_column)
  {
    pch = separate(&buff, par->delim);
    k++;
  }
  while (pch != NULL)
  {
    pi->colnames[j] = pch;
    j++;
    if (j == colnamessize)
    {
      colnamessize *= 2;
      pi->colnames =
          (char **)realloc(pi->colnames, colnamessize * sizeof(char *));
    }
    pch = separate(&buff, par->delim);
  }
  pi->dimensions[1] = j;
  i = 0;
  len = 0;
  nread = get_full_line(&line, &len, datafile, par->gz, gzdatafile);
  while (nread != -1)
  {
    i++;
    free(line);
    line = NULL;
    len = 0;
    nread = get_full_line(&line, &len, datafile, par->gz, gzdatafile);
  }
  Rcpp::Rcout << "done.\n";
  pi->dimensions[0] = i;
  free(line);
  return (0);
}

static int
initialize_dims_h5(struct dim_vars **dvars, struct hdf5_vars **h5vars,
                   struct fastR_user_params *par, struct preprocess_info *pre,
                   struct read_buffers *readbuff, char *placeholder)
{
  size_t bytes, datasize, i, j, remainder;
  float gb, *rowp, **tablep;
  int filecount, file_poi;
  char *name;
  size_t total = 0;
  herr_t status;
  ChunkConfig chunked;
  if (par->transpose == 0)
  {
    // one read buffer per HDF5 file required for transpose==0
    if (par->single == 1)
    {
      file_poi = pre->dimensions[1];
    }
    else if (par->poi_per_file > 0)
    {
      file_poi = par->poi_per_file;
    }
    else
    {
      chunked = Chunker::estimate_num_files(
          (int)pre->dimensions[1], (int)pre->dimensions[0], par->server_threads,
          par->server_memory);
      file_poi = chunked.num_poi;
    }

    if (file_poi < (int)par->chunk_edge)
    {
      par->chunk_edge = file_poi;
    }
    filecount = pre->dimensions[1] / file_poi;
    remainder = pre->dimensions[1] % file_poi;
    if (remainder > 0)
    {
      filecount++;
    }
    Rcpp::Rcout << "POI per H5 file: " << file_poi
                << "\nH5 File Count: " << filecount << "\n";
    pre->row_dim = (par->data_buffer_max - filecount * sizeof(float **)) /
                   (1024 + sizeof(char *) + filecount * sizeof(float *) +
                    pre->dimensions[1] * sizeof(float));
    if (pre->row_dim >= par->chunk_edge)
    {
      pre->row_dim = par->chunk_edge;
    }
    else if (pre->row_dim > 0)
    {
      bytes = filecount * sizeof(float **) +
              par->chunk_edge *
                  (1024 + sizeof(char *) + filecount * sizeof(float *) +
                   pre->dimensions[1] * sizeof(float));
      gb = ((float)bytes) / (1024 * 1024 * 1024);
      gb += 0.01;
      Rprintf("Warning: sub-optimal data buffer size selected. For best "
              "performance specify\n    %.2f Gb or higher\n",
              gb);
    }
    else
    {
      bytes = 1024 + sizeof(char *) + filecount * sizeof(float **) +
              filecount * sizeof(float *) + pre->dimensions[1] * sizeof(float);
      gb = ((float)bytes) / (1024 * 1024 * 1024);
      gb += 0.01;
      Rprintf("Error: selected data buffer size too small to retain a single "
              "row. Select\n    %.2f Gb or higher\n",
              gb);
      return (-1);
    }
    pre->growdim = 0;
    pre->fixdim = 1;
    *dvars = (struct dim_vars *)malloc(filecount * sizeof(struct dim_vars));
    *h5vars = (struct hdf5_vars *)malloc(filecount * sizeof(struct hdf5_vars));
    for (i = 0; (int)i < (filecount - 1); i++)
    {
      (*dvars)[i].vals_memspace_dims[0] = pre->row_dim;
      (*dvars)[i].vals_memspace_dims[1] = file_poi;
      (*dvars)[i].vals_dataspace_dims[0] = pre->dimensions[0];
      (*dvars)[i].vals_dataspace_dims[1] = file_poi;
      (*dvars)[i].vals_chunk_dims[0] = par->chunk_edge;
      (*dvars)[i].vals_chunk_dims[1] = par->chunk_edge;
      (*dvars)[i].vals_hyperslab_pos[0] = 0;
      (*dvars)[i].vals_hyperslab_pos[1] = 0;
    }
    (*dvars)[i].vals_memspace_dims[0] = pre->row_dim;
    (*dvars)[i].vals_dataspace_dims[0] = pre->dimensions[0];
    (*dvars)[i].vals_hyperslab_pos[0] = 0;
    (*dvars)[i].vals_hyperslab_pos[1] = 0;
    if (remainder > 0)
    {
      (*dvars)[i].vals_memspace_dims[1] = remainder;
      (*dvars)[i].vals_dataspace_dims[1] = remainder;
      (*dvars)[i].vals_chunk_dims[0] =
          (par->chunk_edge > remainder ? remainder : par->chunk_edge);
      (*dvars)[i].vals_chunk_dims[1] =
          (par->chunk_edge > remainder ? remainder : par->chunk_edge);
    }
    else
    {
      (*dvars)[i].vals_memspace_dims[1] = file_poi;
      (*dvars)[i].vals_dataspace_dims[1] = file_poi;
      (*dvars)[i].vals_chunk_dims[0] = par->chunk_edge;
      (*dvars)[i].vals_chunk_dims[1] = par->chunk_edge;
    }
    // allocate value and row name read buffers
    datasize = filecount * sizeof(float **) +
               filecount * pre->row_dim * sizeof(float *) +
               pre->row_dim * pre->dimensions[1] * sizeof(float);
    readbuff->val_buffer = (float ***)malloc(datasize);
    tablep = (float **)(readbuff->val_buffer + filecount);
    rowp = (float *)(tablep + filecount * pre->row_dim);
    for (i = 0; (int)i < filecount; i++)
    {
      readbuff->val_buffer[i] = (float **)(tablep + i * pre->row_dim);
      for (j = 0; j < pre->row_dim; j++)
      {
        readbuff->val_buffer[i][j] =
            (float *)(rowp + i * pre->row_dim * file_poi +
                      j * (*dvars)[i].vals_memspace_dims[1]);
      }
    }
    readbuff->row_buffer = (char **)malloc(pre->row_dim * sizeof(char *));
    for (j = 0; j < pre->row_dim; j++)
    {
      readbuff->row_buffer[j] = placeholder;
    }
  }
  else
  {
    // only one read buffer required when transpose==1
    if (par->single == 1)
    {
      file_poi = pre->dimensions[0];
    }
    else if (par->poi_per_file > 0)
    {
      file_poi = par->poi_per_file;
    }
    else
    {
      chunked = Chunker::estimate_num_files(
          (int)pre->dimensions[0], (int)pre->dimensions[1], par->server_threads,
          par->server_memory);
      file_poi = chunked.num_poi;
    }
    if (file_poi < (int)par->chunk_edge)
    {
      par->chunk_edge = file_poi;
    }
    filecount = pre->dimensions[0] / file_poi;
    remainder = pre->dimensions[0] % file_poi;
    if (remainder > 0)
    {
      filecount++;
    }
    Rcpp::Rcout << "POI per H5 file: " << file_poi
                << "\nH5 File Count: " << filecount << "\n";
    pre->row_dim = (par->data_buffer_max - sizeof(float **)) /
                   (1024 + sizeof(char *) + sizeof(float *) +
                    pre->dimensions[1] * sizeof(float));
    if (pre->row_dim >= par->chunk_edge)
    {
      pre->row_dim = par->chunk_edge;
    }
    else if (pre->row_dim > 0)
    {
      bytes = sizeof(float **) +
              par->chunk_edge * (1024 + sizeof(char *) + sizeof(float *) +
                                 pre->dimensions[1] * sizeof(float));
      gb = ((float)bytes) / (1024 * 1024 * 1024);
      gb += 0.01;
      Rprintf("Warning: sub-optimal data buffer size selected. For best "
              "performance specify\n    %.2f Gb or higher\n",
              gb);
    }
    else
    {
      bytes = sizeof(float **) + 1024 + sizeof(char *) + sizeof(float *) +
              pre->dimensions[1] * sizeof(float);
      gb = ((float)bytes) / (1024 * 1024 * 1024);
      gb += 0.01;
      Rprintf("Error: selected data buffer size too small to retain a single "
              "row. Select\n    %.2f Gb or higher\n",
              gb);
      return (-1);
    }
    pre->growdim = 1;
    pre->fixdim = 0;
    *dvars = (struct dim_vars *)malloc(filecount * sizeof(struct dim_vars));
    *h5vars = (struct hdf5_vars *)malloc(filecount * sizeof(struct hdf5_vars));
    // flip row and column dimension if matrix is transposed
    for (i = 0; (int)i < (filecount - 1); i++)
    {
      (*dvars)[i].vals_memspace_dims[0] = pre->dimensions[1];
      (*dvars)[i].vals_memspace_dims[1] = pre->row_dim;
      (*dvars)[i].vals_dataspace_dims[0] = pre->dimensions[1];
      (*dvars)[i].vals_dataspace_dims[1] = file_poi;
      (*dvars)[i].vals_chunk_dims[0] = par->chunk_edge;
      (*dvars)[i].vals_chunk_dims[1] = par->chunk_edge;
      (*dvars)[i].vals_hyperslab_pos[0] = 0;
      (*dvars)[i].vals_hyperslab_pos[1] = 0;
    }
    (*dvars)[i].vals_memspace_dims[0] = pre->dimensions[1];
    (*dvars)[i].vals_dataspace_dims[0] = pre->dimensions[1];
    (*dvars)[i].vals_hyperslab_pos[0] = 0;
    (*dvars)[i].vals_hyperslab_pos[1] = 0;
    if (remainder > 0)
    {
      (*dvars)[i].vals_dataspace_dims[1] = remainder;
      (*dvars)[i].vals_memspace_dims[1] = remainder;
      (*dvars)[i].vals_chunk_dims[0] =
          (par->chunk_edge > remainder ? remainder : par->chunk_edge);
      (*dvars)[i].vals_chunk_dims[1] =
          (par->chunk_edge > remainder ? remainder : par->chunk_edge);
    }
    else
    {
      (*dvars)[i].vals_dataspace_dims[1] = file_poi;
      (*dvars)[i].vals_memspace_dims[1] = pre->row_dim;
      (*dvars)[i].vals_chunk_dims[0] = par->chunk_edge;
      (*dvars)[i].vals_chunk_dims[1] = par->chunk_edge;
    }
    // allocate value and row name read buffers
    datasize = sizeof(float **) + pre->dimensions[1] * sizeof(float *) +
               pre->dimensions[1] * pre->row_dim * sizeof(float);
    readbuff->val_buffer = (float ***)malloc(datasize);
    tablep = (float **)(readbuff->val_buffer + 1);
    rowp = (float *)(tablep + pre->dimensions[1]);
    readbuff->val_buffer[0] = tablep;
    for (j = 0; j < pre->dimensions[1]; j++)
    {
      readbuff->val_buffer[0][j] = (float *)(rowp + j * pre->row_dim);
    }
    readbuff->row_buffer = (char **)malloc(pre->row_dim * sizeof(char *));
    for (j = 0; j < pre->row_dim; j++)
    {
      readbuff->row_buffer[j] = placeholder;
    }
  }
  // create output HDF5 files
  Rcpp::Rcout << "Initializing H5 files..." << std::endl;
  for (i = 0; (int)i < filecount; i++)
  {
    // create file
    name = (char *)malloc(2 * strlen(par->h5file_base) + 100);
    Rcpp::Rcout << name << par->h5file_base << "/" << par->h5file_base << "." << i << ".h5" << std::endl;
    (*h5vars)[i].file =
        H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    free(name);
    // create row dataspace
    (*h5vars)[i].row_dataspace = H5Screate_simple(
        1, &((*dvars)[i].vals_dataspace_dims[pre->growdim]), NULL);
    (*h5vars)[i].row_memspace = H5Screate_simple(
        1, &((*dvars)[i].vals_memspace_dims[pre->growdim]), NULL);
    (*h5vars)[i].row_datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size((*h5vars)[i].row_datatype, H5T_VARIABLE);
    (*h5vars)[i].row_prop_list = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk((*h5vars)[i].row_prop_list, 1,
                 &((*dvars)[i].vals_chunk_dims[0]));
    H5Pset_deflate((*h5vars)[i].row_prop_list, 4);
    (*h5vars)[i].row_dataset =
        H5Dcreate((*h5vars)[i].file, pre->rowlabel, (*h5vars)[i].row_datatype,
                  (*h5vars)[i].row_dataspace, H5P_DEFAULT,
                  (*h5vars)[i].row_prop_list, H5P_DEFAULT);
    // create column dataset and write to all files
    (*h5vars)[i].col_dataspace = H5Screate_simple(
        1, &((*dvars)[i].vals_dataspace_dims[pre->fixdim]), NULL);
    (*h5vars)[i].col_datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size((*h5vars)[i].col_datatype, H5T_VARIABLE);
    (*h5vars)[i].col_prop_list = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk((*h5vars)[i].col_prop_list, 1,
                 &((*dvars)[i].vals_chunk_dims[0]));
    H5Pset_deflate((*h5vars)[i].col_prop_list, 4);
    (*h5vars)[i].col_dataset =
        H5Dcreate((*h5vars)[i].file, pre->collabel, (*h5vars)[i].col_datatype,
                  (*h5vars)[i].col_dataspace, H5P_DEFAULT,
                  (*h5vars)[i].col_prop_list, H5P_DEFAULT);
    status =
        H5Dwrite((*h5vars)[i].col_dataset, (*h5vars)[i].col_datatype, H5S_ALL,
                 H5S_ALL, H5P_DEFAULT, &(pre->colnames[total * pre->fixdim]));
    total += (*dvars)[i].vals_dataspace_dims[pre->fixdim];
    if (status < 0)
    {
      Rprintf("Error: unable to write column names to HDF5 file\n");
      return (-1);
    }
    // create values dataset
    (*h5vars)[i].vals_prop_list = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk((*h5vars)[i].vals_prop_list, 2, (*dvars)[i].vals_chunk_dims);
    H5Pset_deflate((*h5vars)[i].vals_prop_list, 4);
    (*h5vars)[i].vals_dataspace =
        H5Screate_simple(2, (*dvars)[i].vals_dataspace_dims, NULL);
    (*h5vars)[i].vals_memspace =
        H5Screate_simple(2, (*dvars)[i].vals_memspace_dims, NULL);
    (*h5vars)[i].vals_datatype = H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_order((*h5vars)[i].vals_datatype, H5T_ORDER_LE);
    (*h5vars)[i].vals_dataset =
        H5Dcreate((*h5vars)[i].file, "values", (*h5vars)[i].vals_datatype,
                  (*h5vars)[i].vals_dataspace, H5P_DEFAULT,
                  (*h5vars)[i].vals_prop_list, H5P_DEFAULT);
    (*h5vars)[i].written = 0;
  }
  Rcpp::Rcout << "done.\n";
  return (filecount);
}

static int read_write_rownames_values(
    FILE *datafile, gzFile gzdatafile, struct read_buffers *readbuff,
    struct hdf5_vars *h5vars, struct dim_vars *dv,
    struct fastR_user_params *par, int filecount, struct preprocess_info *pre)
{
  size_t f, fmin, x, y, i, j, k, l, iadj, idim, jdim, vcfbuffmax, *indptr,
      *predptr, *fptr;
  float **dstarts, val;
  hsize_t position_zero[2] = {0, 0};
  ssize_t nread;
  char *line, *pch = NULL;
  char *buff, *vcfbuff;
  size_t len = 0;
  size_t skip = 0;
  size_t fdummy = 0;
  herr_t status;
  int done = 0;
  Rcpp::Rcout << "Writing H5 files..." << std::endl;
  // flip data buffer index variable if matrix is transposed
  if (par->transpose == 1)
  {
    indptr = &j;
    predptr = &i;
    fptr = &fdummy;
    idim = 1;
    jdim = 0;
  }
  else
  {
    indptr = &i;
    predptr = &j;
    fptr = &f;
    idim = 0;
    jdim = 1;
  }
  // H5Dwrite needs pointer to start of data matrix rather than block of row
  // pointers
  dstarts = (float **)malloc(filecount * sizeof(float *));
  for (f = 0; (int)f < filecount; f++)
  {
    dstarts[f] = &readbuff->val_buffer[*fptr][0][0];
  }
  // read input file line by line and fill data buffer
  f = 0;
  fmin = 0;
  x = 0;
  y = 0;
  i = 0;
  iadj = par->header_row;
  for (k = 0; k <= par->header_row; k++)
  {
    nread = get_full_line(&line, &len, datafile, par->gz, gzdatafile);
    if (nread == -1)
    {
      Rprintf("Error: could not read first data row from input file\n");
      free(dstarts);
      free(line);
      return (1);
    }
  }
  buff = line;
  while (nread != -1)
  {
    if (done != 0)
    {
      Rprintf("Warning: additional data found at end of input file,\n  "
                      " contents may have changed during processing\n");
      free(dstarts);
      free(line);
      return (0);
    }
    j = 0;
    y = 0;
    k = 0;
    f = fmin;
    if (par->vcf == 0)
    {
      while (k < par->name_column)
      {
        pch = separate(&buff, par->delim);
        k++;
      }
      readbuff->row_buffer[i] = strdup(pch);
      while (k < par->data_column)
      {
        pch = separate(&buff, par->delim);
        k++;
      }
    }
    else
    {
      vcfbuffmax = 100;
      vcfbuff = (char *)malloc(vcfbuffmax * sizeof(char));
      *vcfbuff = '\0';
      while (k < par->data_column)
      {
        pch = separate(&buff, par->delim);
        k++;
        switch (k)
        {
        case 1:
        case 2:
        case 4:
        case 5:
          while (vcfbuffmax < (strlen(vcfbuff) + strlen(pch) + 2))
          {
            vcfbuffmax = 2 * vcfbuffmax;
            vcfbuff = (char *)realloc(vcfbuff, vcfbuffmax * sizeof(char));
          }
          strcat(vcfbuff, pch);
          if (k < 5)
          {
            strcat(vcfbuff, "_");
          }
          else
          {
            // check for presense of multiple ALTs, if so, skip to next line
            l = 0;
            while (*(pch + l) != '\0')
            {
              if (*(pch + l) == ',')
              {
                skip++;
                free(vcfbuff);
                vcfbuff = NULL;
                f++;
                break;
              }
              l++;
            }
          }
          break;
        }
      }
      if (vcfbuff == NULL)
      {
        free(line);
        line = NULL;
        len = 0;
        nread = get_full_line(&line, &len, datafile, par->gz, gzdatafile);
        buff = line;
        continue;
      }
      readbuff->row_buffer[i] = strdup(vcfbuff);
      free(vcfbuff);
    }
    while (pch != NULL && y < pre->dimensions[1])
    {
      // check explicitly for empty field at the end of each line
      if (pch[0] == '\0' || pch[0] == '\n')
      {
        readbuff->val_buffer[*fptr][*indptr][*predptr] = NAN;
      }
      else
      {
        if (par->vcf == 0)
        {
          val = (float)atof(pch);
          // atof returns 0 if no valid float is identified - in that case,
          // ensure that the value begins with a numeral, otherwise write NAN to
          // buffer
          if (val == 0.0 &&
              (pch[0] < 43 || pch[0] > 57 || pch[0] == 47 || pch[0] == 44))
          {
            val = NAN;
          }
        }
        else
        {
          // genotype is missing if not 0/0, 0|0, 0|1, 1|0, 0/1, 1|1, or 1/1
          switch ((char)*(pch + 1))
          {
          case '|':
          case '/':
            switch ((int)*pch + (int)*(pch + 2) - 96)
            {
            case 0:
              val = 0.0;
              break;
            case 1:
              val = 1.0;
              break;
            case 2:
              val = 2.0;
              break;
            default:
              val = NAN;
              break;
            }
            break;
          default:
            if (*(pch + 1) == '\0')
            {
              switch ((int)*pch - 48)
              {
              case 0:
                val = 0.0;
                break;
              case 1:
                val = 1.0;
                break;
              default:
                val = NAN;
                break;
              }
            }
            else
            {
              val = NAN;
            }
            break;
          }
        }
        readbuff->val_buffer[*fptr][*indptr][*predptr] = val;
      }
      pch = separate(&buff, par->delim);
      j++;
      y++;
      if (j == dv[f].vals_memspace_dims[jdim])
      {
        j = 0;
        f++;
      }
    }
    i++;
    x++;
    // ensure each line contains the expect field count; fill remainder of
    // buffer with NAN
    if (j != 0)
    {
      Rprintf("Warning: fields missing in input file line %zu\n",
              iadj + i + skip);
      i -= 1;
      while (j < dv[f].vals_memspace_dims[jdim])
      {
        readbuff->val_buffer[*fptr][*indptr][*predptr] = NAN;
        j++;
      }
      j = 0;
      f++;
      i++;
    }
    else if (pch != NULL)
    {
      Rprintf(
          "Warning: more fields than expected found in input file line %zu\n",
          iadj + i + skip);
    }
    // write to HDF5 file once buffer is full or current file/files are full
    if (i == pre->row_dim || x == dv[fmin].vals_dataspace_dims[idim])
    {
      for (k = fmin; k < f; k++)
      {
        dv[k].vals_memspace_dims[idim] = i;
        h5vars[k].written = 1;
        H5Sselect_hyperslab(h5vars[k].vals_memspace, H5S_SELECT_SET,
                            position_zero, NULL, dv[k].vals_memspace_dims,
                            NULL);
        H5Sselect_hyperslab(h5vars[k].vals_dataspace, H5S_SELECT_SET,
                            dv[k].vals_hyperslab_pos, NULL,
                            dv[k].vals_memspace_dims, NULL);
        status = H5Dwrite(h5vars[k].vals_dataset, H5T_NATIVE_FLOAT,
                          h5vars[k].vals_memspace, h5vars[k].vals_dataspace,
                          H5P_DEFAULT, dstarts[k]);
        if (status < 0)
        {
          Rprintf("Error: unable to write data block to HDF5 file\n");
          free(dstarts);
          free(line);
          return (1);
        }
        H5Sselect_hyperslab(h5vars[k].row_memspace, H5S_SELECT_SET,
                            &(position_zero[0]), NULL,
                            &(dv[k].vals_memspace_dims[pre->growdim]), NULL);
        H5Sselect_hyperslab(h5vars[k].row_dataspace, H5S_SELECT_SET,
                            &(dv[k].vals_hyperslab_pos[pre->growdim]), NULL,
                            &(dv[k].vals_memspace_dims[pre->growdim]), NULL);
        status = H5Dwrite(h5vars[k].row_dataset, h5vars[k].row_datatype,
                          h5vars[k].row_memspace, h5vars[k].row_dataspace,
                          H5P_DEFAULT, readbuff->row_buffer);
        if (status < 0)
        {
          Rprintf("Error: unable to write row name block to HDF5 file\n");
          free(dstarts);
          free(line);
          return (1);
        }
        dv[k].vals_hyperslab_pos[pre->growdim] += i;
      }
      // clean up row name buffer before refilling - sizes of identifiers may
      // vary
      for (k = 0; k < i; k++)
      {
        free(readbuff->row_buffer[k]);
        readbuff->row_buffer[k] = NULL;
      }
      iadj += i;
      i = 0;
      if (x == dv[fmin].vals_dataspace_dims[idim])
      {
        x = 0;
        fmin++;
        if ((int)par->transpose == 0 || ((int)par->transpose == 1 && (int)fmin == filecount))
        {
          done = 1;
        }
      }
    }
    free(line);
    line = NULL;
    len = 0;
    nread = get_full_line(&line, &len, datafile, par->gz, gzdatafile);
    buff = line;
  }
  // at end of file, check for any unwritten data, write and modify data set
  // dimensions as necessary
  if (done == 0)
  {
    l = 0;
    for (k = fmin; (int)k < filecount; k++)
    {
      l += dv[k].vals_dataspace_dims[idim];
    }
    if (l != (x + skip))
    {
      Rcpp::Rcout << "Warning: fewer data entries found in input file than expected, \n"
                  << "contents may have changed during processing" << std::endl;
    }
    if (i > 0)
    {
      for (k = fmin; k < f; k++)
      {
        h5vars[k].written = 1;
        dv[k].vals_memspace_dims[idim] = i;
        dv[k].vals_dataspace_dims[idim] = x;
        H5Dset_extent(h5vars[k].vals_dataset, dv[k].vals_dataspace_dims);
        H5Dset_extent(h5vars[k].row_dataset,
                      &dv[k].vals_dataspace_dims[pre->growdim]);
        H5Sselect_hyperslab(h5vars[k].vals_memspace, H5S_SELECT_SET,
                            position_zero, NULL, dv[k].vals_memspace_dims,
                            NULL);
        H5Sselect_hyperslab(h5vars[k].vals_dataspace, H5S_SELECT_SET,
                            dv[k].vals_hyperslab_pos, NULL,
                            dv[k].vals_memspace_dims, NULL);
        status = H5Dwrite(h5vars[k].vals_dataset, H5T_NATIVE_FLOAT,
                          h5vars[k].vals_memspace, h5vars[k].vals_dataspace,
                          H5P_DEFAULT, dstarts[k]);
        if (status < 0)
        {
          Rprintf("Error: unable to write data block to HDF5 file\n");
          free(dstarts);
          free(line);
          return (1);
        }
        H5Sselect_hyperslab(h5vars[k].row_memspace, H5S_SELECT_SET,
                            &(position_zero[0]), NULL,
                            &(dv[k].vals_memspace_dims[pre->growdim]), NULL);
        H5Sselect_hyperslab(h5vars[k].row_dataspace, H5S_SELECT_SET,
                            &dv[k].vals_hyperslab_pos[pre->growdim], NULL,
                            &dv[k].vals_memspace_dims[pre->growdim], NULL);
        status = H5Dwrite(h5vars[k].row_dataset, h5vars[k].row_datatype,
                          h5vars[k].row_memspace, h5vars[k].row_dataspace,
                          H5P_DEFAULT, readbuff->row_buffer);
        if (status < 0)
        {
          Rprintf("Error: unable to write row name block to HDF5 file\n");
          free(dstarts);
          free(line);
          return (1);
        }
        dv[k].vals_hyperslab_pos[pre->growdim] += i;
      }
    }
  }
  Rcpp::Rcout << "done." << std::endl;
  free(line);
  free(dstarts);
  return (0);
}

static void preprocess_cleanup(struct preprocess_info *pre)
{
  free(pre->collabel);
  free(pre->rowlabel);
  free(pre->colnames);
  free(pre->colline);
  return;
}

static void final_cleanup(struct read_buffers *readbuff,
                          struct hdf5_vars **h5vars, struct dim_vars **dv,
                          struct preprocess_info *pi,
                          struct fastR_user_params *par, int filecount)
{
  size_t i;
  char *name;
  for (i = 0; i < pi->row_dim; i++)
  {
    free(readbuff->row_buffer[i]);
  }
  free(readbuff->row_buffer);
  free(readbuff->val_buffer);
  for (i = 0; (int)i < filecount; i++)
  {
    H5Sclose((*h5vars)[i].vals_dataspace);
    H5Sclose((*h5vars)[i].col_dataspace);
    H5Sclose((*h5vars)[i].row_dataspace);
    H5Sclose((*h5vars)[i].row_memspace);
    H5Sclose((*h5vars)[i].vals_memspace);
    H5Tclose((*h5vars)[i].vals_datatype);
    H5Tclose((*h5vars)[i].col_datatype);
    H5Tclose((*h5vars)[i].row_datatype);
    H5Dclose((*h5vars)[i].vals_dataset);
    H5Dclose((*h5vars)[i].col_dataset);
    H5Dclose((*h5vars)[i].row_dataset);
    H5Fclose((*h5vars)[i].file);
    if ((*h5vars)[i].written == 0)
    {
      name = (char *)malloc(2 * strlen(par->h5file_base) + 100);
      Rcpp::Rcout << name << par->h5file_base << "/" << par->h5file_base << "." << i << ".h5" << std::endl;
      // sprintf(name, "%s/%s.%03d.h5", par->h5file_base, par->h5file_base,
      //         (int)i);
      remove(name);
      free(name);
    }
  }
  free(*h5vars);
  free(*dv);
  return;
}

static int execute_fastR_hdf5convert(struct fastR_user_params *up)
{
  FILE *infile = NULL;
  gzFile gzinfile = NULL;
  struct hdf5_vars *h5v;
  struct dim_vars *dims;
  struct read_buffers rb;
  struct preprocess_info pre;
  int flag;
  int nfiles;
  char *dummy;
  dummy = strdup("X");
  if (up->vcf == 1)
  {
    up->transpose = 1;
    up->data_column = 10;
    free(up->delim);
    up->delim = strdup("\t");
  }
  if (up->gz == 1)
  {
    gzinfile = gzopen(up->infile_name, "r");
    if (gzinfile == NULL)
    {
      Rprintf("Error: cannot open input file \"%s\"\n",
              up->infile_name);
      return (1);
    }
    flag = preprocess_datafile(infile, gzinfile, up, &pre);
    if (flag == 1)
    {
      preprocess_cleanup(&pre);
      return (1);
    }
    gzclose(gzinfile);
    gzinfile = gzopen(up->infile_name, "r");
  }
  else
  {
    infile = fopen(up->infile_name, "r");
    if (infile == NULL)
    {
      Rprintf("Error: cannot open input file \"%s\"\n",
              up->infile_name);
      return (1);
    }
    flag = preprocess_datafile(infile, gzinfile, up, &pre);
    if (flag == 1)
    {
      preprocess_cleanup(&pre);
      return (1);
    }
    fclose(infile);
    infile = fopen(up->infile_name, "r");
  }
  nfiles = initialize_dims_h5(&dims, &h5v, up, &pre, &rb, dummy);
  preprocess_cleanup(&pre);
  free(dummy);
  if (nfiles == -1)
  {
    final_cleanup(&rb, &h5v, &dims, &pre, up, nfiles);
    return (1);
  }
  flag = read_write_rownames_values(infile, gzinfile, &rb, h5v, dims, up,
                                    nfiles, &pre);
  final_cleanup(&rb, &h5v, &dims, &pre, up, nfiles);
  if (up->gz == 1)
  {
    gzclose(gzinfile);
  }
  else
  {
    fclose(infile);
  }
  if (flag == 1)
  {
    return (1);
  }
  return (0);
}

// [[Rcpp::export]]
int FastRegImportCpp(std::string dataFile, std::string h5File, int headerRow,
                     int idCol, int dataCol, float buffSize,
                     bool transpose, int chunkEdge, bool vcf,
                     std::string delim, bool gz, int poiPerFile,
                     bool singleFile, int serverThreads, float serverMem)
{
  struct fastR_user_params par;
  int retval;
  par.infile_name = strdup(dataFile.c_str());
  par.h5file_base = strdup(h5File.c_str());
  par.poi_per_file = poiPerFile;
  par.header_row = (size_t)headerRow;
  par.name_column = (size_t)idCol;
  par.data_column = (size_t)dataCol;
  par.data_buffer_max = (size_t)(buffSize * 1024 * 1024 * 1024);
  par.transpose = (int)transpose;
  par.chunk_edge = (hsize_t)chunkEdge;
  par.vcf = (int)vcf;
  par.delim = strdup(delim.c_str());
  par.gz = (int)gz;
  par.single = (int)singleFile;
  par.server_threads = serverThreads;
  par.server_memory = serverMem;
  retval = execute_fastR_hdf5convert(&par);
  free(par.infile_name);
  free(par.h5file_base);
  free(par.delim);
  return (retval);
}
