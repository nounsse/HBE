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
 * @file util.h
 * @brief util header 
 * @author Yi-Qing WANG <yiqing.wang@polytechnique.edu>
 */


#ifndef EIGEN_H
#define EIGEN_H
#include <Eigen/Dense>  	// Eigen
#endif
using namespace Eigen;

/** @brief coordinate of a pixel in an image */
struct SamplingParams{
	int patch_size;
	double flat_threshold;
	double orient_threshold;
	int num_orientations;
	int flat_model;
	int textural_model;
	int minimal_num_samples;
	bool verbose;
	const char * filename;
};

/** @brief exception messaging */
void fail (const char * message );

/** @brief routine timer */
void times (const char * which );

/**
 * @brief read in a PNG and return its dimension 
 *
 * @param file_name the image name
 * @param image the image as a MatrixXd array 
 */
void imread(	
		const char * file_name
,		int & image_rows
,		int & image_cols
, 		int & num_channels
,		MatrixXd *& image		
);

void imread(size_t nx, size_t ny, size_t nc,
        float * pixel_stream,
        MatrixXd *& image);

void imread_16bits(	
		const char * file_name
,		int & image_rows
,		int & image_cols
, 		int & num_channels
,		MatrixXd *& image		
);

void imread_channs_16b(	
		const char * file_name
,		int i_num_channels
,		int & image_rows
,		int & image_cols
, 		int & num_channels
,		MatrixXd *& image		
);


/**
 * @brief write out as a PNG
 *
 * @param file_name the image name
 * @param num_channels RGB or Gray
 */
void imwrite( 
		const char * file_name 
, 		MatrixXd * image
, 		int num_channels
);

void imwrite_pgm( 
		const char * file_name 
, 		MatrixXd * image
, 		int num_channels
);

/**
 * @brief reduce a gray image into patches 
 *
 * @param image a matrix representing one image
 * @param overlap number of columns (resp. rows) two horizontal (resp. vectical) neighbor patches share
 * @param patch_size patch dimension: patch_size * patch_size 
 * @return an array of vectorized patches 
 */
VectorXd ** image2patches( 
		MatrixXd const * image
,		int image_rows
,		int image_cols
,		int num_channels
,		int overlap
,		int patch_size 
);


/** @brief number of patches required to cover the whole image by row or by column */
int num_patches( 
		int n_pixels
,		int patch_size
,		int overlap
);


/**
 * @brief the inverse of image2patches 
 *
 * @param patch_at_coordinates an array of vectorized patches 
 * @param image the image made of the patches in the patch array 
 * @param normalize whether to average the estimates on each pixel 
 */
MatrixXd * patches2image(     
		VectorXd ** patch_at_coordinates
,		int overlap
, 		int patch_size
,		int num_channels
, 		int image_rows
, 		int image_cols
//,		double * likelihoods_winners
,		bool normalize = true
);


/**
 * @brief switch between RBG and a slightly different YUV 
 *
 * @param image image array. In case of a gray image, do nothing 
 * @param inverse RGB to YUV or YUV to RGB 
 */
void RGB_transform( 
		MatrixXd * image
,		int num_channels
,		bool inverse 
);


/**@brief present the patch array with matrices*/
void vectorArray2Matrix( 
		VectorXd ** varray
,		int row_id
,		MatrixXd & vmatrix
);


void separate_channels( float *u0, float **u1, int ncol, int nrow );

MatrixXd * remake_bayer( MatrixXd * image, int num_channels, int num_channels_out );
