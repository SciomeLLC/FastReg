#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>

//consolidates all HDF5 handles and attributes
struct hdf5_vars {
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
};

//consolidates all dimension variables
struct dim_vars {
	hsize_t vals_dataspace_dims[2];
	hsize_t vals_chunk_dims[2];
	hsize_t vals_memspace_dims[2];
	hsize_t vals_hyperslab_pos[2];
	hsize_t unlimited_dims[2];
	hsize_t unlimited_single_dim;
	hsize_t col_dim;
	hsize_t row_dim;
	hsize_t row_chunk_dim;
	int growdim;
};

//consolidates user-defined parameters
struct fastR_user_params {
	char *infile_name;
	char *h5file_name;
	size_t header_row;
	size_t name_column;
	size_t data_column;
	size_t data_buffer_max;
	int transpose;
	hsize_t chunk_edge;
};

//consolidates read buffer pointers
struct read_buffers {
	float **val_buffer;
	char **row_buffer;
};

//implements POSIX function strsep for portability
char *separate(char **restrict stringp, const char *restrict delim) {
	size_t span;
	char *first, *last;
	if(*stringp==NULL) {
		return(NULL);
	}
	first=*stringp;
	span=strcspn(first, delim);
	last=first+span;
	if(*last=='\0') {
		*stringp=NULL;
	} else {
		*last='\0';
		*stringp=last+1;
	}
	return(first);
}

//implements POSIX function getline for portability
ssize_t get_full_line(char **restrict lineptr, size_t *restrict n, FILE *restrict stream) {
	char *res, *start, *mark;
	size_t i, prev=0, count;
	if(*lineptr==NULL || *n==0) {
		*lineptr=(char*)malloc(10485760);
		*n=10485760;
		if(*lineptr==NULL) {
			return(-1);
		}
	}
	start=*lineptr;
	count=*n;
	res=fgets(start, (int)count, stream);
	if(res==NULL) {
		return(-1);
	} else {
		while(1) {
			mark=(char*)memchr((void*)start, '\n', count);
			if(mark!=NULL) {
				return((ssize_t)(mark-*lineptr)+1);
			}
			prev=*n;
			*n=2*prev;
			*lineptr=(char*)realloc(*lineptr, *n);
			if(*lineptr==NULL) {
				return(-1);
			}
			start=*lineptr+prev-1;
			count=prev+1;
			res=fgets(start, (int)count, stream);
			if(res==NULL) {
				return(-1);
			}
		}
	}
}

//initialize all dimensions
static void init_dims(struct dim_vars *dv, int chunk) {
	dv->vals_dataspace_dims[0]=0;
	dv->vals_dataspace_dims[1]=0;
	dv->vals_chunk_dims[0]=chunk;
	dv->vals_chunk_dims[1]=chunk;
	dv->vals_memspace_dims[0]=0;
	dv->vals_memspace_dims[1]=0;
	dv->vals_hyperslab_pos[0]=0;
	dv->vals_hyperslab_pos[1]=0;
	dv->unlimited_dims[0]=H5S_UNLIMITED;
	dv->unlimited_dims[1]=H5S_UNLIMITED;
	dv->unlimited_single_dim=H5S_UNLIMITED;
	dv->col_dim=0;
	dv->row_dim=0;
	dv->row_chunk_dim=chunk;
	dv->growdim=0;
	return;
}

static int read_header_write_cols(FILE *datafile, struct hdf5_vars *h5vars, struct dim_vars *dv, struct fastR_user_params *par) {
	ssize_t nread;
	char *line=NULL;
	char *buff, *pch=NULL;
	char **colnames, *collabel;
	size_t colnamessize;
	size_t len=0;
	hsize_t i, j, k;
	herr_t status;
	i=0;
	j=0;
	k=0;
	//set column data set name based on matrix orientation
	if(par->transpose==1) {
		collabel=strdup("individuals");
	} else {
		collabel=strdup("predictors_of_interest");
	}
	colnames=(char**)malloc(sizeof(char*));
	colnamessize=1;
	//toss lines until hear row is read
	nread=get_full_line(&line, &len, datafile);
	i++;
	if(nread==-1) {
		fprintf(stderr, "Error: could not read specified header row from input file\n");
		return(1);
	}
	while(i<par->header_row) {
		free(line);
		line=NULL;
		len=0;
		nread=get_full_line(&line, &len, datafile);
		i++;
		if(nread==-1) {
			fprintf(stderr, "Error: could not read specified header row from input file\n");
			return(1);
		}
	}
	//parse and store column names
	if(line[nread-1]=='\n') {
		line[nread-1]='\0';
	}
	//hold line pointer so memory can be freed after parsing
	buff=line;
	//toss fields before first data point
	while(k<par->data_column) {
		pch=separate(&buff,"\t ");
		k++;
	}
	while(pch!=NULL) {
		colnames[j]=pch;
		j++;
		if(j==colnamessize) {
			colnamessize*=2;
			colnames=(char**)realloc(colnames, colnamessize*sizeof(char*));
		}
		pch=separate(&buff,"\t ");
	}
	dv->col_dim=j;
	//write column names to HDF5 file
	h5vars->col_dataspace=H5Screate_simple(1, &dv->col_dim, &dv->unlimited_single_dim);
	h5vars->col_datatype=H5Tcopy(H5T_C_S1);
	H5Tset_size(h5vars->col_datatype, H5T_VARIABLE);
	h5vars->col_prop_list=H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(h5vars->col_prop_list, 1, &dv->row_chunk_dim);
	H5Pset_deflate(h5vars->col_prop_list, 4);
	h5vars->col_dataset=H5Dcreate(h5vars->file, collabel, h5vars->col_datatype, h5vars->col_dataspace, H5P_DEFAULT, h5vars->col_prop_list, H5P_DEFAULT);
	status=H5Dwrite(h5vars->col_dataset, h5vars->col_datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, colnames);
	if(status<0) {
		fprintf(stderr, "Error: unable to write column names to HDF5 file\n");
		return(1);
	}
	free(collabel);
	free(colnames);
	free(line);
	return(0);
}

static int allocate_buffers_finalize_dims(struct read_buffers *readbuff, struct dim_vars *dv, struct fastR_user_params *par) {
	//finalize dimensions
	size_t bytes, datasize, i;
	float gb, *rowp;
	dv->row_dim=par->data_buffer_max/(1024+sizeof(char*)+sizeof(float*)+dv->col_dim*sizeof(float));
	if(dv->row_dim>=100) {
		dv->row_dim=100;
	} else if(dv->row_dim>0) {
		bytes=100*(1024+sizeof(char*)+sizeof(float*)+dv->col_dim*sizeof(float));
		gb=((float)bytes)/(1024*1024*1024);
		gb+=0.01;
		fprintf(stderr, "Warning: sub-optimal data buffer size selected. For best performance specify\n    %.2f Gb or higher\n", gb);
	} else {
		bytes=1024+sizeof(char*)+sizeof(float*)+dv->col_dim*sizeof(float);
		gb=((float)bytes)/(1024*1024*1024);
		gb+=0.01;
		fprintf(stderr, "Error: selected data buffer size too small to retain a single row. Select\n    %.2f Gb or higher\n", gb);
		return(1);
	}
	//flip row and column dimension if matrix is transposed
	if(par->transpose==1) {
		dv->vals_dataspace_dims[0]=dv->col_dim;
		dv->vals_dataspace_dims[1]=dv->row_dim;
		dv->vals_memspace_dims[0]=dv->col_dim;
		dv->vals_memspace_dims[1]=dv->row_dim;
		dv->growdim=1;
	} else {
		dv->vals_dataspace_dims[0]=dv->row_dim;
		dv->vals_dataspace_dims[1]=dv->col_dim;
		dv->vals_memspace_dims[0]=dv->row_dim;
		dv->vals_memspace_dims[1]=dv->col_dim;
		dv->growdim=0;
	}
	//allocate full read buffers
	datasize=dv->vals_memspace_dims[0]*sizeof(float*)+dv->vals_memspace_dims[0]*dv->vals_memspace_dims[1]*sizeof(float);
	readbuff->val_buffer=(float**)malloc(datasize);
	rowp=(float*)(readbuff->val_buffer+dv->vals_memspace_dims[0]);
	for(i=0;i<dv->vals_memspace_dims[0];i++) {
		readbuff->val_buffer[i]=rowp+dv->vals_memspace_dims[1]*i;
	}
	readbuff->row_buffer=(char**)malloc(dv->row_dim*sizeof(char*));
	return(0);
}

static void create_rownames_values_datasets(struct hdf5_vars *h5vars, struct dim_vars *dv, struct fastR_user_params *par) {
	char *rowlabel;
	//initialize row name dataset
	if(par->transpose==1) {
		rowlabel=strdup("predictors_of_interest");
	} else {
		rowlabel=strdup("individuals");
	}
	h5vars->row_dataspace=H5Screate_simple(1, &dv->row_dim, &dv->unlimited_single_dim);
	h5vars->row_memspace=H5Screate_simple(1, &dv->row_dim, NULL);
	h5vars->row_datatype=H5Tcopy(H5T_C_S1);
	H5Tset_size(h5vars->row_datatype, H5T_VARIABLE);
	h5vars->row_prop_list=H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(h5vars->row_prop_list, 1, &dv->row_chunk_dim);
	H5Pset_deflate(h5vars->row_prop_list, 4);
	h5vars->row_dataset=H5Dcreate(h5vars->file, rowlabel, h5vars->row_datatype, h5vars->row_dataspace, H5P_DEFAULT, h5vars->row_prop_list, H5P_DEFAULT);
	free(rowlabel);
	//initialize values data set
	h5vars->vals_prop_list=H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(h5vars->vals_prop_list, 2, dv->vals_chunk_dims);
	H5Pset_deflate(h5vars->vals_prop_list, 4);
	h5vars->vals_dataspace=H5Screate_simple(2, dv->vals_dataspace_dims, dv->unlimited_dims);
	h5vars->vals_memspace=H5Screate_simple(2, dv->vals_memspace_dims, NULL);
	h5vars->vals_datatype=H5Tcopy(H5T_NATIVE_FLOAT);
	H5Tset_order(h5vars->vals_datatype, H5T_ORDER_LE);
	h5vars->vals_dataset=H5Dcreate(h5vars->file, "values", h5vars->vals_datatype, h5vars->vals_dataspace, H5P_DEFAULT, h5vars->vals_prop_list, H5P_DEFAULT);
	return;
}

static int read_write_rownames_values(FILE *datafile, struct read_buffers *readbuff, struct hdf5_vars *h5vars, struct dim_vars *dv, struct fastR_user_params *par) {
	size_t i, j, k, iadj, *indptr, *predptr;
	float *dstart, val;
	ssize_t nread;
	char *line=NULL;
	char *buff, *pch;
	size_t len=0;
	herr_t status;
	//flip data buffer index variable if matrix is transposed
	if(par->transpose==1) {
		indptr=&j;
		predptr=&i;
	} else {
		indptr=&i;
		predptr=&j;
	}
	//H5Dwrite needs pointer to start of data matrix rather than block of row pointers
	dstart=&readbuff->val_buffer[0][0];
	//read input file line by line and fill data buffer
	i=0;
	nread=get_full_line(&line, &len, datafile);
	if(nread==-1) {
		fprintf(stderr, "Error: could not read first data row from input file\n");
		return(1);
	}
	buff=line;
	while(nread!=-1) {
		j=0;
		k=0;
		while(k<par->name_column) {
			pch=separate(&buff,"\t ");
			k++;
		}
		readbuff->row_buffer[i]=strdup(pch);
		while(k<par->data_column) {
			pch=separate(&buff,"\t ");
			k++;
		}
		//check explicitly for empty field at the end of each line
		while(pch!=NULL) {
			if(pch[0]=='\0' || pch[0]=='\n') {
				readbuff->val_buffer[*indptr][*predptr]=NAN;
			}
			else {
				val=(float)atof(pch);
				//atof returns 0 if no valid float is identified - in that case, ensure that the value begins with a numeral, otherwise write NAN to buffer
				if(val==0.0 && (pch[0]<43 || pch[0]>57 || pch[0]==47 || pch[0]==44)) {
					val=NAN;
				}
				readbuff->val_buffer[*indptr][*predptr]=val;
			}
			pch=separate(&buff,"\t ");
			j++;
		}
		i++;
		//ensure each line contains the expect field count
		if(j!=dv->col_dim) {
			iadj=i+par->header_row;
			fprintf(stderr, "Error: fields missing in input file line %zu\n", iadj);
		}
		//write to HDF5 file once buffer is full
		if(i==dv->row_dim) {
			H5Sselect_hyperslab(h5vars->vals_dataspace, H5S_SELECT_SET, dv->vals_hyperslab_pos, NULL, dv->vals_memspace_dims, NULL);
			status=H5Dwrite(h5vars->vals_dataset, H5T_NATIVE_FLOAT, h5vars->vals_memspace, h5vars->vals_dataspace, H5P_DEFAULT, dstart);
			if(status<0) {
				fprintf(stderr, "Error: unable to write data block to HDF5 file\n");
				return(1);
			}
			H5Sselect_hyperslab(h5vars->row_dataspace, H5S_SELECT_SET, &dv->vals_hyperslab_pos[dv->growdim], NULL, &dv->vals_memspace_dims[dv->growdim], NULL);
			status=H5Dwrite(h5vars->row_dataset, h5vars->row_datatype, h5vars->row_memspace, h5vars->row_dataspace, H5P_DEFAULT, readbuff->row_buffer);
			if(status<0) {
				fprintf(stderr, "Error: unable to write row name block to HDF5 file\n");
				return(1);
			}
			dv->vals_hyperslab_pos[dv->growdim]+=dv->row_dim;
			dv->vals_dataspace_dims[dv->growdim]+=dv->row_dim;
			H5Dset_extent(h5vars->vals_dataset, dv->vals_dataspace_dims);
			H5Sset_extent_simple(h5vars->vals_dataspace, 2, dv->vals_dataspace_dims, dv->unlimited_dims);
			H5Dset_extent(h5vars->row_dataset, &dv->vals_dataspace_dims[dv->growdim]);
			H5Sset_extent_simple(h5vars->row_dataspace, 1, &dv->vals_dataspace_dims[dv->growdim], &dv->unlimited_single_dim);
			//clean up row name buffer before refilling - sizes of identifiers may vary
			for(k=0;k<dv->vals_memspace_dims[dv->growdim];k++) {
				free(readbuff->row_buffer[k]);
			}
			i=0;
		}
		free(line);
		line=NULL;
		len=0;
		nread=get_full_line(&line, &len, datafile);
		buff=line;
	}
	//when end of file is reached, shrink dimensions to match final block to be written
	dv->vals_memspace_dims[dv->growdim]=i%dv->row_dim;
	dv->vals_dataspace_dims[dv->growdim]=dv->vals_dataspace_dims[dv->growdim]-dv->row_dim+dv->vals_memspace_dims[dv->growdim];
	H5Dset_extent(h5vars->vals_dataset, dv->vals_dataspace_dims);
	H5Sset_extent_simple(h5vars->vals_dataspace, 2, dv->vals_dataspace_dims, dv->vals_dataspace_dims);
	H5Sset_extent_simple(h5vars->vals_memspace, 2, dv->vals_memspace_dims, dv->vals_memspace_dims);
	H5Dset_extent(h5vars->row_dataset, &dv->vals_dataspace_dims[dv->growdim]);
	H5Sset_extent_simple(h5vars->row_dataspace, 1, &dv->vals_dataspace_dims[dv->growdim], &dv->vals_dataspace_dims[dv->growdim]);
	H5Sset_extent_simple(h5vars->row_memspace, 1, &dv->vals_memspace_dims[dv->growdim], &dv->vals_memspace_dims[dv->growdim]);
	H5Sselect_hyperslab(h5vars->vals_dataspace, H5S_SELECT_SET, dv->vals_hyperslab_pos, NULL, dv->vals_memspace_dims, NULL);
	status=H5Dwrite(h5vars->vals_dataset, H5T_NATIVE_FLOAT, h5vars->vals_memspace, h5vars->vals_dataspace, H5P_DEFAULT, dstart);
	if(status<0) {
		fprintf(stderr, "Error: unable to write final data block to HDF5 file\n");
		return(1);
	}
	H5Sselect_hyperslab(h5vars->row_dataspace, H5S_SELECT_SET, &dv->vals_hyperslab_pos[dv->growdim], NULL, &dv->vals_memspace_dims[dv->growdim], NULL);
	status=H5Dwrite(h5vars->row_dataset, h5vars->row_datatype, h5vars->row_memspace, h5vars->row_dataspace, H5P_DEFAULT, readbuff->row_buffer);
	if(status<0) {
		fprintf(stderr, "Error: unable to write final row name block to HDF5 file\n");
		return(1);
	}
	free(line);
	return(0);
}

static void final_cleanup(struct read_buffers *readbuff, struct hdf5_vars *h5vars, struct dim_vars *dv) {
	size_t i;
	for(i=0;i<dv->vals_memspace_dims[dv->growdim];i++) {
		free(readbuff->row_buffer[i]);
	}
	free(readbuff->row_buffer);
	H5Sclose(h5vars->vals_dataspace);
	H5Sclose(h5vars->col_dataspace);
	H5Sclose(h5vars->row_dataspace);
	H5Sclose(h5vars->row_memspace);
	H5Sclose(h5vars->vals_memspace);
	H5Tclose(h5vars->vals_datatype);
	H5Tclose(h5vars->col_datatype);
	H5Tclose(h5vars->row_datatype);
	H5Dclose(h5vars->vals_dataset);
	H5Dclose(h5vars->col_dataset);
	H5Dclose(h5vars->row_dataset);
	free(readbuff->val_buffer);
	H5Fclose(h5vars->file);
	return;
}

static int execute_fastR_hdf5convert(struct fastR_user_params *up) {
	FILE *infile;
	struct hdf5_vars h5v;
	struct dim_vars dims;
	struct read_buffers rb;
	int flag;
	init_dims(&dims, up->chunk_edge);
	infile=fopen(up->infile_name, "r");
	if(infile==NULL) {
		fprintf(stderr, "Error: cannot open input file \"%s\"\n", up->infile_name);
		return(1);
	}
	h5v.file=H5Fcreate(up->h5file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	flag=read_header_write_cols(infile, &h5v, &dims, up);
	if(flag==1) {
		return(1);
	}
	flag=allocate_buffers_finalize_dims(&rb, &dims, up);
	if(flag==1) {
		return(1);
	}
	create_rownames_values_datasets(&h5v, &dims, up);
	flag=read_write_rownames_values(infile, &rb, &h5v, &dims, up);
	if(flag==1) {
		return(1);
	}
	final_cleanup(&rb, &h5v, &dims);
	fclose(infile);
	return(0);
}

SEXP fastR_hdf5convert(SEXP dataFile, SEXP h5File, SEXP headerRow, SEXP idCol, SEXP dataCol, SEXP buffSize, SEXP transpose, SEXP chunkEdge) {
	struct fastR_user_params par;
	int retval;
	SEXP result=PROTECT(allocVector(LGLSXP, 1));
	par.infile_name=strdup(CHAR(asChar(dataFile)));
	par.h5file_name=strdup(CHAR(asChar(h5File)));
	strcat(par.h5file_name, ".h5");
	par.header_row=asInteger(headerRow);
	par.name_column=asInteger(idCol);
	par.data_column=asInteger(dataCol);
	par.data_buffer_max=(size_t)(asReal(buffSize)*1024*1024*1024);
	par.transpose=asInteger(transpose);
	par.chunk_edge=asInteger(chunkEdge);
	retval=execute_fastR_hdf5convert(&par);
	if(retval==0) {
		LOGICAL(result)[0]=1;
	} else {
		LOGICAL(result)[0]=0;
	}
	free(par.infile_name);
	free(par.h5file_name);
	UNPROTECT(1);
	return(result);
}
