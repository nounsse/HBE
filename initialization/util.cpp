/*
 * Copyright 2009-2012 Yi-Qing WANG 
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/**
 * @file util.cpp
 * @brief utils for both PLE and SPLE 
 * @author Yi-Qing WANG <yiqing.wang@polytechnique.edu>
 */

#include <sys/time.h>		// time the routines
#include <fstream>		// read in noise data
#include "io_png.h" 		// image IO
#include "io_pgm.h" 		// image IO
#include "io_exr.h"
#include <stack>

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef UTIL_H
#define UTIL_H
#include "util.h"
#endif


#ifndef RAND_H
#define RAND_H
extern "C"{
	#include "randmt.h"	//external random number generator
}
#endif

#ifndef STD_H
#define STD_H
#include <iomanip>
#include <iostream>
#include <cstdio>
#endif

#ifndef DATA_H
#define DATA_H
#include "DataProvider.h" 	// initial mixture
#endif

#ifndef SHORT_NEWMAT
#define SHORT_NEWMAT
#include "../newmat10/newmatap.h"  // NEWMAT
typedef NEWMAT::Matrix NMatrix;
typedef NEWMAT::SymmetricMatrix NSym;
typedef NEWMAT::DiagonalMatrix NDiag;
typedef NEWMAT::ColumnVector NColV;
typedef NEWMAT::RowVector NRowV;
#endif

typedef Matrix<double,Dynamic,Dynamic,RowMajor> RMatrixXd;
using namespace std;

void fail (const char *message) {
	cerr << message << endl;
	exit(EXIT_FAILURE);
}

//routine timer
void times (const char *which) {
	/* If which is not empty, print the times since the previous call. */
	static double last_wall = 0.0, last_cpu = 0.0;
	double wall, cpu;
	struct timeval tv;
	clock_t stamp;
	
	wall = last_wall;
	cpu = last_cpu;
	if (gettimeofday(&tv,NULL) != 0 || (stamp = clock()) == (clock_t)-1)
	    fail("Unable to get times");
	last_wall = tv.tv_sec+1.0e-6*tv.tv_usec;
	last_cpu = stamp/(double)CLOCKS_PER_SEC;
	if (strlen(which) > 0){
	    wall = last_wall-wall;
	    cpu = last_cpu-cpu;
	    printf("%s time = %.2f seconds, CPU = %.2f seconds\n",which,wall,cpu);
	}
}


//wrapper in Matlab style
void imread(	
		const char * file_name
,		int & image_rows
,		int & image_cols
, 		int & num_channels
,		MatrixXd *& image		
){
		size_t nx, ny, nc;
		float * pixel_stream = NULL;
		pixel_stream = read_png_f32( file_name, &nx, &ny, &nc );
		if( pixel_stream == NULL )
			fail("Unable to get the image");	
		
		//return these useful parameters 
		image_cols = (int)nx;
		image_rows = (int)ny;
		num_channels = (int)nc;
		
		//input stream assumes row-major while Eigen defaults to column-major
		Map<MatrixXf> parallel( pixel_stream, image_cols, image_rows * num_channels );
		image = new MatrixXd [ num_channels ];
		for( int ch = 0; ch < num_channels; ch++ )
			image[ ch ] = parallel.block( 0, ch*image_rows, image_cols, image_rows).transpose().cast<double>();

		//release
		free(pixel_stream);
		cout << "INFO: read in image " << file_name << " of dimension: (" << ny << " , " << nx << " , " << nc << ")" << endl;
}

//wrapper in Matlab style
void imread_16bits(	
		const char * file_name
,		int & image_rows
,		int & image_cols
, 		int & num_channels
,		MatrixXd *& image		
){
		int nx, ny, nc;
		float * pixel_stream = NULL;
		nc = 1;
		pixel_stream = read_pgm_float(file_name, &nx, &ny);	
		if( pixel_stream == NULL )
			fail("Unable to get the image");	
		
		//return these useful parameters 
		image_cols = (int)nx;
		image_rows = (int)ny;
		num_channels = (int)nc;
		
		//input stream assumes row-major while Eigen defaults to column-major
		Map<MatrixXf> parallel( pixel_stream, image_cols, image_rows * num_channels );
		image = new MatrixXd [ num_channels ];
		for( int ch = 0; ch < num_channels; ch++ )
			image[ ch ] = parallel.block( 0, ch*image_rows, image_cols, image_rows).transpose().cast<double>();

		//release
		free(pixel_stream);
		cout << "INFO: read in image " << file_name << " of dimension: (" << ny << " , " << nx << " , " << nc << ")" << endl;
}

//wrapper in Matlab style
void imread_channs_16b(	
		const char * file_name
,		int i_num_channels
,		int & image_rows
,		int & image_cols
, 		int & num_channels
,		MatrixXd *& image		
){
		int nx, ny, nc;
		float * pixel_stream = NULL;
		nc = i_num_channels;
		pixel_stream = read_pgm_float(file_name, &nx, &ny);	
		if( pixel_stream == NULL )
			fail("Unable to get the image");			
		
		if ( i_num_channels == 4 ) {
			image_cols = (int)nx / 2;
			image_rows = (int)ny / 2;
		} else if ( i_num_channels == 1 ) {
			image_cols = (int)nx;
			image_rows = (int)ny;
		} else 			
			fail("Incorrect number of channels. Code done for 1 or 4 channels. Use 1 for gray level images and 4 for RAW images (color images).");			
		
		num_channels = (int)nc;
		
		float *fpI[num_channels];
		for (int ii=0; ii < num_channels; ii++) 
			fpI[ii] = (float*) malloc(image_rows*image_cols*sizeof(float)); 
		
		if ( i_num_channels == 4 )
			separate_channels( pixel_stream, fpI, nx, ny );
		else if ( i_num_channels == 1 )
			for ( int k=0; k < image_rows*image_cols; k++)
				fpI[0][k] = pixel_stream[k];
		else 
			cout << "Incorrect number of channels\n" << endl;

		//input stream assumes row-major while Eigen defaults to column-major
		//Map<MatrixXf> parallel( pixel_stream, image_cols, image_rows * num_channels );
		image = new MatrixXd [ num_channels ];
		for( int ch = 0; ch < num_channels; ch++ ) {
			Map<MatrixXf> parallel( fpI[ch], image_cols, image_rows );
			image[ ch ] = parallel.transpose().cast<double>();
		}
		//release
		free(pixel_stream);
		cout << "INFO: read in image " << file_name << " of dimension: (" << ny << " , " << nx << " , " << nc << ")" << endl;
}

//write image and release memory
void imwrite( 
		const char * file_name 
, 		MatrixXd * image
, 		int num_channels
){
	int image_cols = image[0].cols();
	int image_rows = image[0].rows();
	int pixels_per_channel = image_cols * image_rows; 
	float * output = new float [ pixels_per_channel * num_channels ];
	//this part should be straightforward but still be careful with the order
	#pragma omp parallel for schedule( static )
	for(int j = 0; j < image_cols; j++)
		for(int i = 0; i < image_rows; i++)
			for(int ch = 0; ch < num_channels; ch++)
				output[ ch*pixels_per_channel + i*image_cols + j ] = (float) image[ ch ](i,j);
	//release
	delete [] image;
	write_png_f32( file_name, output, (size_t) image_cols, (size_t) image_rows, (size_t) num_channels );
	delete [] output;
	cout << "INFO: write the image " << file_name << " to local folder." << endl;
}

//write image and release memory
void imwrite_pgm( 
		const char * file_name 
, 		MatrixXd * image
, 		int num_channels
){
	int image_cols = image[0].cols();
	int image_rows = image[0].rows();
	int pixels_per_channel = image_cols * image_rows; 
	float * output = new float [ pixels_per_channel * num_channels ];
	//this part should be straightforward but still be careful with the order
	#pragma omp parallel for schedule( static )
	for(int j = 0; j < image_cols; j++)
		for(int i = 0; i < image_rows; i++)
			for(int ch = 0; ch < num_channels; ch++)
				output[ ch*pixels_per_channel + i*image_cols + j ] = (float) image[ ch ](i,j);
	//release
	delete [] image;
	write_pgm_float( file_name, output, image_cols, image_rows );
	delete [] output;
	cout << "INFO: write the image " << file_name << " to local folder." << endl;
}

//write image and release memory
void imwrite_exr( 
		const char * file_name 
, 		MatrixXd * image
, 		int num_channels
){
	int image_cols = image[0].cols();
	int image_rows = image[0].rows();
	int pixels_per_channel = image_cols * image_rows; 
	
	float * output = new float [ pixels_per_channel * num_channels ];

	//this part should be straightforward but still be careful with the order
	#pragma omp parallel for schedule( static )
	for(int j = 0; j < image_cols; j++)
		for(int i = 0; i < image_rows; i++)
			for(int ch = 0; ch < num_channels; ch++)
				output[ ch*pixels_per_channel + i*image_cols + j ] = (float) image[ch](i,j);
	//release
	//delete [] image;	
	writeEXR_float( file_name, output, image_cols, image_rows );
	delete [] output;
	cout << "INFO: write the image " << file_name << " to local folder." << endl;
}

void imread_exr(	
		const char * file_name
,		int i_num_channels
,		int & image_rows
,		int & image_cols
, 		int & num_channels
,		MatrixXd *& image		
){
		int nx, ny, nc;
		float * pixel_stream = NULL;
		nc = i_num_channels;
		pixel_stream = readEXR_float(file_name, &nx, &ny);	
		if( pixel_stream == NULL )
			fail("Unable to get the image.");			
		
		if ( i_num_channels == 4 ) {
			image_cols = (int)nx / 2;
			image_rows = (int)ny / 2;
		} else if ( i_num_channels == 1 ) {
			image_cols = (int)nx;
			image_rows = (int)ny;
		} else 			
			fail("Incorrect number of channels. Code done for 1 or 4 channels. Use 1 for gray level images and 4 for RAW images (color images).");			
		
		num_channels = (int)nc;
		
		float *fpI[num_channels];
		for (int ii=0; ii < num_channels; ii++) 
			fpI[ii] = (float*) malloc(image_rows*image_cols*sizeof(float)); 
		
		if ( i_num_channels == 4 )
			separate_channels( pixel_stream, fpI, nx, ny );
		else if ( i_num_channels == 1 )
			for ( int k=0; k < image_rows*image_cols; k++)
				fpI[0][k] = pixel_stream[k];
		else 
			cout << "Incorrect number of channels\n" << endl;

		//input stream assumes row-major while Eigen defaults to column-major
		//Map<MatrixXf> parallel( pixel_stream, image_cols, image_rows * num_channels );
		image = new MatrixXd [ num_channels ];
		for( int ch = 0; ch < num_channels; ch++ ) {
			Map<MatrixXf> parallel( fpI[ch], image_cols, image_rows );
			image[ ch ] = parallel.transpose().cast<double>();
		}
		//release
		free(pixel_stream);
		cout << "INFO: read in image " << file_name << " of dimension: (" << ny << " , " << nx << " , " << nc << ")" << endl;
}


//add noise to image
void add_gaussian_noise( 
		MatrixXd * image
, 		double sigma
,		int num_channels
,		bool generate_noise 
){
	int image_rows = image[0].rows();
	int image_cols = image[0].cols();
	int pixels_per_channel = image_rows * image_cols;

	//generate defaults to true
	//the other option allows you to read in noise
	//which is helpful in algo comparison and debugging
	if( generate_noise ){
		init_randmt_auto();
		double * noise = new double [ pixels_per_channel ];
		for(int ch = 0; ch < num_channels; ch++){
			for( int i = 0; i < pixels_per_channel; i++ )
				noise[ i ] = sigma * rand_normal();
			Map<MatrixXd> pure_noise( noise, image_rows, image_cols );
			image[ch] += pure_noise;
		}
		delete [] noise;
	}else{
		//noise_data is a binary file containing noise
		int num_data = pixels_per_channel * num_channels;
		double noise [num_data];
		FILE * noise_data = fopen( "noise_data", "rb" );
		if( noise_data == NULL )
			fail("Unable to get the noise data");
		size_t count = fread( noise, sizeof(double), num_data, noise_data );
		Map<MatrixXd> pure_noise( noise, image_rows, image_cols * num_channels );
		for( int ch = 0; ch < num_channels; ch++ )
			image[ch] += pure_noise.block( 0, ch * image_cols, image_rows, image_cols );
	}
	
	cout << "INFO: noise with sigma = " << sigma << " is added to the image." << endl;
	
}

//calc num of patches needed to cover the whole image
//according to a sliding window scheme specified by overlap
//which is the num of columns shared by horizontally neighboring patches
int num_patches( 
		int n_pixels
,		int patch_size
,		int overlap
){
	int step = patch_size - overlap;
//	it holds that for some k, n_pixels = patch_size + step * k + something
//	with something = 0 to k-1
	int something = (n_pixels - patch_size) % step;
	int correction;
	if ( something == 0 )
		correction = 1;
	else
		correction = 2;
	int k = (n_pixels - something - patch_size)/step;
	return k + correction;
}


//reduce an image into patches
VectorXd ** image2patches( 
		MatrixXd const * image
,		int image_rows
,		int image_cols
,		int num_channels
,		int overlap
,		int patch_size 
){
	int map_rows = num_patches( image_rows, patch_size, overlap );
	int map_cols = num_patches( image_cols, patch_size, overlap );
	int data_size = pow( patch_size, 2 );
	/*cout << "in function " << endl;
	cout << image_rows << " " << image_cols << " " <<  patch_size << " " <<  overlap << endl; 
	cout << "Vals: " << map_rows << " " << map_cols << endl;*/
	cout << "overlap img2pt " << overlap << endl;
	
//	allocate some memory for patches
	VectorXd ** patch_at_coordinates = new VectorXd * [ map_rows ];
	for( int row = 0; row < map_rows; row++ )
		patch_at_coordinates[ row ] = new VectorXd [ map_cols ];
//	the patch upper left corner's coordinates
	int coordinate_j, coordinate_i = -1*patch_size;
//	coordinate cannot exceed max 
	int max_coordinate_i = image_rows - patch_size;
	int max_coordinate_j = image_cols - patch_size;
//	sliding window step
	int step = patch_size - overlap;
//	meat
	
	for( int i = 0; i < map_rows; i++ ){
		coordinate_i = max( 0, min( max_coordinate_i, coordinate_i + step ) );
		coordinate_j = -1*patch_size;		
		for( int j = 0; j < map_cols; j++ ){
			patch_at_coordinates[i][j] = VectorXd::Zero(data_size*num_channels);
			coordinate_j = max( 0, min( max_coordinate_j, coordinate_j + step ) );
			for ( int ch=0; ch < num_channels; ch++ ) {
				/*for (unsigned kk=0; kk<patch_size; kk++)
					for (unsigned jj=0; jj<patch_size; jj++)
						cout << "img = " << image[ch]( coordinate_i+kk, coordinate_j+jj) << " ";
				cout<< endl;*/
				MatrixXd patch (patch_size, patch_size);
				patch = image[ch].block( coordinate_i, coordinate_j, patch_size, patch_size ).transpose();								
				patch_at_coordinates[i][j].segment(ch*data_size,data_size) = Map<MatrixXd>( patch.data(), data_size, 1 );
			}
		}
	}
	cout << "INFO: the image is reduced to ( " << map_rows <<  " , " << map_cols << " ) patches." << endl;
	return patch_at_coordinates;
}

//assemble patches to form an image again
MatrixXd * patches2image(     
		VectorXd ** patch_at_coordinates
,		int overlap
, 		int patch_size
,		int num_channels
, 		int image_rows
, 		int image_cols
//,		double *likelihoods_winners
,		bool normalize
){
	int data_size = patch_size * patch_size;
	MatrixXd * image = new MatrixXd [num_channels];
//  reset the whole picture
	for (int ch=0; ch < num_channels; ch++)
	    image[ch].setZero( image_rows, image_cols );
//  mask counts for each pixel the number of patches covering it
	MatrixXd * mask = new MatrixXd [num_channels];
//  reset the whole picture
	for (int ch=0; ch < num_channels; ch++)
	    mask[ch].setZero( image_rows, image_cols );
    int coordinate_j, coordinate_i = -1*patch_size;
	int max_coordinate_i = image_rows - patch_size;
	int max_coordinate_j = image_cols - patch_size;
	int step = patch_size - overlap;
//	block to mark the patch in the mask 
	MatrixXd block( patch_size, patch_size );
	block.setOnes( patch_size, patch_size );
    int map_rows = num_patches( image_rows, patch_size, overlap );
    int map_cols = num_patches( image_cols, patch_size, overlap );

    for( int i = 0; i < map_rows; i++ ){
                coordinate_i = max( 0, min( max_coordinate_i,  coordinate_i + step ) );
                coordinate_j = -1*patch_size;
                for( int j = 0; j < map_cols; j++ ){
                        coordinate_j = max( 0, min( max_coordinate_j, coordinate_j + step ) );
                        for ( int ch=0; ch < num_channels; ch++ ) {
			//a transposition to reflect what has been done in image2patches 
							VectorXd patch_aux = patch_at_coordinates[i][j].segment(ch*data_size,data_size);							
	                        image[ch].block(coordinate_i, coordinate_j, patch_size, patch_size) += Map<MatrixXd>( patch_aux.data(), patch_size, patch_size ).transpose();// * likelihoods_winners[ i*map_cols + j];
	                        
			//equivalently, use RowMajor config 
                        //image.block(coordinate_i, coordinate_j, patch_size, patch_size) += Map<RMatrixXd>( patch_at_coordinates[i][j].data(), patch_size, patch_size ); 
			if( normalize )
	                        mask[ch].block(coordinate_i, coordinate_j, patch_size, patch_size) += block; //*likelihoods_winners[ i*map_cols + j];
						}
                }
    }
	if( normalize )
		for ( int ch=0; ch < num_channels; ch++ ) 
			image[ch] = image[ch].cwiseQuotient( mask[ch] );

	return image;
}


void vectorArray2Matrix( 
		VectorXd ** varray
,		int row_id
,		MatrixXd & vmatrix
){
	int n_cols = vmatrix.cols();
	int n_rows = vmatrix.rows();
	for( int col = 0; col < n_cols; col++ )
		vmatrix.block( 0, col, n_rows, 1 ) = varray[row_id][col];
}

//tensor structure orientation detector
int tensorStructure(
		MatrixXd & image
,		int row
,		int col
,		SamplingParams & params	
){
	Matrix2d tensor = Matrix2d::Zero();
	static int patch_size = params.patch_size;
	//first I need to compute gradient at all the pixel sites in the patch
	for( int i = row; i < row + patch_size; i++ )	
		for( int j = col; j < col + patch_size; j++ ){
			//I'll always avoid image boundary, so relax
			//this scheme of differentiation is better as the second order error is gone
			double row_grad = image(i+1,j)	- image(i-1,j);
			double col_grad = image(i,j+1) - image(i,j-1);
			Vector2d grad;
			grad << row_grad, col_grad;	
			tensor += grad*grad.transpose();
		}
	//the traditional Ax = \lambda x thing didn't prove numerically reliable. Sorry, fall back to NewMat
	NMatrix cov_container(2, 2);
	cov_container << tensor.data(); 
	NSym sym_cov(2);
	sym_cov << cov_container;
	NMatrix U(2, 2);
	NDiag D(2);
	SVD( sym_cov, D, U  );
	double big_val = D(1,1);
	double small_val = D(2,2); 
	//now we are going to decide 
	int model;
	if( big_val/small_val < params.orient_threshold ){
		model = big_val < params.flat_threshold ? params.flat_model : params.textural_model;
	}else{
		//take the eigenvector associated with the larger eigenvalue
		double theta = atan(U(2,1)/U(1,1));
		theta = theta >= 0 ? theta : theta + M_PI;
		model = floor(theta*params.num_orientations/M_PI);
		//sometimes, if theta is a very small negative, shit does happen
		if( model == params.textural_model ){
			cout << "INFO: a tiny numerical approximation leads to a big quantification error. Corrected!" << endl;
			model = model - 1;
		}	
	}
	return model;
}


void separate_channels( float *u0, float **u1, int ncol, int nrow ) 
{
  /* We suppose that the number of rows and columns of the input images 
   * are even. Otherwise, the size of the output images would be different
   * for the different channels. */
	
  if ( nrow % 2 != 0 ) /* number of rows is not even */
    printf("Warning: the number of rows has been truncated to be even.\n");

  if ( ncol % 2 != 0 ) /* number of columns is not even */
    printf("Warning: the number of columns has been truncated to be even.\n");
	
  int nrow_ch = nrow / 2;	
  int ncol_ch = ncol / 2;	
  int k,x, y = 0;
	
  for ( k=0; k < nrow_ch; k++, y++ ) {
    int l = y*ncol;
    int h = k*ncol_ch;
		
    float *ptr = &u0[l];
    float *ptrR = &u1[0][h];
    float *ptrG1 = &u1[1][h];
		
    for ( x=0; x < ncol_ch; x++, ptr++, ptrR++, ptrG1++ ) {
      *ptrR = *ptr;   /* Fill red channel */				
      ptr++;
      *ptrG1 = *ptr;  /* Fill green channel */				
    }		

    y++;
    l = y*ncol;
    ptr = &u0[l];
    float *ptrG2 = &u1[2][h];
    float *ptrB = &u1[3][h];
		
    for ( x=0; x < ncol_ch; x++, ptr++, ptrG2++, ptrB++ ) {
      *ptrG2 = *ptr;   /* Fill green channel */
      ptr++;
      *ptrB = *ptr;    /* Fill blue channel */
    }

  }
	
}

MatrixXd * remake_bayer( MatrixXd * image, int num_channels, int num_channels_out ) 
{
	
  /* u1 is the input 4 channel image */
  /* u0 is the reconstructed Bayer pattern */
	
  /* We suppose that the number of rows and columns of the input images 
   * are even. Otherwise, the size of the output images would be different
   * for the different channels. */
  int ncol_ch = image[0].cols();
  int nrow_ch = image[0].rows();
  
  int ncol = ncol_ch * 2;	
  int wh = ncol_ch*nrow_ch;
  int k,x,y = 0;
  
  float* u1[num_channels];
  for ( int ch=0; ch<num_channels; ch++ ) {
	u1[ch] = (float*)malloc(wh*sizeof(float));
	for ( int i=0; i<nrow_ch ; i++ )
		for ( int j=0; j< ncol_ch; j++ ) {
			//MatrixXd imageT = image[ch].transpose();
			u1[ch][i*ncol_ch + j] = (float)image[ch](i,j);
		}
  }			
  
  float* u0 = (float *)malloc(4*ncol_ch*nrow_ch*sizeof(float));
	
  for ( k=0; k < nrow_ch; k++, y++ ) {
    int l = y*ncol;
    int h = k*ncol_ch;
		
    float *ptr = &u0[l];
    float *ptrR = &u1[0][h];
    float *ptrG1 = &u1[1][h];
		
    for ( x=0; x < ncol_ch; x++, ptr++, ptrR++, ptrG1++ ) {
      *ptr = *ptrR;   /* Fill red channel */
      ptr++;
      *ptr = *ptrG1;  /* Fill green channel */
    }

    y++;
    l = y*ncol;
    ptr = &u0[l];
    float *ptrG2 = &u1[2][h];
    float *ptrB = &u1[3][h];
		
    for ( x=0; x < ncol_ch; x++, ptr++, ptrG2++, ptrB++ ) {
      *ptr = *ptrG2;  /* Fill green channel */
      ptr++;
      *ptr = *ptrB;   /* Fill blue channel */
    }

  }
  
  MatrixXd * out_image = new MatrixXd [num_channels_out];
  for ( int ch=0; ch<num_channels_out; ch++)	{
	MatrixXd tmp(2*nrow_ch,2*ncol_ch);
	out_image[ch] = tmp;
  }
  for ( int i=0; i<2*nrow_ch; i++)
		for ( int j=0; j<2*ncol_ch; j++)
			out_image[0](i,j) = (double)u0[i*2*ncol_ch+j];
  //out_image[0] = out_image[0].transpose();			
  	
  return out_image;	
}

std::string getFileExt(const std::string& s)
{
    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos)
    {
        return(s.substr(i+1, s.length() - i));
    }
    else
        return("");
}

