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
 * @file PLE_lib.cpp
 * @brief library for PLE denoising algorithm
 * @author Yi-Qing WANG <yiqing.wang@polytechnique.edu>
 */

#ifndef PLE_H
#define PLE_H
#include "PLE_lib.h"
#endif

#include <stack>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef SHORT_NEWMAT
#define SHORT_NEWMAT
#include "../newmat10/newmatap.h"  // NEWMAT
typedef NEWMAT::Matrix NMatrix;
typedef NEWMAT::SymmetricMatrix NSym;
typedef NEWMAT::DiagonalMatrix NDiag;
typedef NEWMAT::ColumnVector NColV;
typedef NEWMAT::RowVector NRowV;
typedef NEWMAT::LogAndSign NLogAndSign;
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
#include <iostream>
#include <iomanip>
#include <cstdio>
#endif
using namespace Eigen;
using namespace std; 

//main_routine adds noise and stores the clean image for comparison purpose
int PLE_main_routine( 	
        const char * input
        ,			const char * output
        ,			const char * u_mtx
        ,			const char * f_mtx
        ,			const char * prnu_mtx
        ,			int overlap
        ,			int patch_size
        ,			int num_orientations
        ,			double muR
        ,			double g
        ,			double sigma2R
        ,			double epsilon
        ,			int i_num_channels
        ){
    int num_channels, image_rows, image_cols;
    MatrixXd * clean_image = NULL;
    imread_channs_16b( input, i_num_channels, image_rows, image_cols, num_channels, clean_image );
    if( clean_image == NULL )
        fail("Unable to get the image");

    int Unum_channels, Uimage_rows, Uimage_cols;
    MatrixXd * U = NULL;
    if (strcmp((const char*)(getFileExt(u_mtx).c_str()),"png")==0){
        imread( u_mtx, Uimage_rows, Uimage_cols,
                Unum_channels, U );
    } else {
        std::cout << "16b " << std::endl;
        imread_channs_16b( u_mtx, i_num_channels, Uimage_rows, Uimage_cols, Unum_channels, U );
    }
    if( U == NULL )
        fail("Unable to get the U image");
    if( Unum_channels != num_channels )
        fail("Number of channels for the mask should match the one of the input");

    int Fnum_channels, Fimage_rows, Fimage_cols;
    MatrixXd * F = NULL;
    imread_exr( f_mtx, i_num_channels, Fimage_rows, Fimage_cols, Fnum_channels, F );
    if( F == NULL )
        fail("Unable to get the F image");

    int Pnum_channels, Pimage_rows, Pimage_cols;
    MatrixXd * PRNU = NULL;
    imread_exr( prnu_mtx, i_num_channels, Pimage_rows, Pimage_cols, Pnum_channels, PRNU );
    if( PRNU == NULL )
        fail("Unable to get the PRNU image");

    MatrixXd * noisy = new MatrixXd [num_channels];
    MatrixXd * invVar = new MatrixXd [num_channels];
    for ( int ch = 0; ch < num_channels; ch++ ) {
        MatrixXd tmp (image_rows,image_cols);
        MatrixXd tmp2 (image_rows,image_cols);
        noisy[ch] = tmp;
        invVar[ch] = tmp2;
        for ( int i=0; i<image_rows; i++)
            for ( int j=0; j<image_cols; j++) {
                noisy[ch](i,j) = U[ch](i,j)*(clean_image[ch](i,j) - muR)/(F[ch](i,j)*PRNU[ch](i,j)*g);
                if (!isfinite(noisy[ch](i,j))) {
                    cout <<  U[ch](i,j) << " " << clean_image[ch](i,j) << " " << muR << " " << F[ch](i,j) << " " << PRNU[ch](i,j) << " " << g << endl;
                    cout << "Warning nan value in noisy values!!!! " << endl;
                    exit (EXIT_FAILURE);
                }
                invVar[ch](i,j) = U[ch](i,j)*(F[ch](i,j)*PRNU[ch](i,j)*g)*(F[ch](i,j)*PRNU[ch](i,j)*g)/(g*(clean_image[ch](i,j) - muR) + sigma2R) + (1-U[ch](i,j))/sigma2R;
            }
    }

    PLE_sub_routine( noisy, invVar, num_channels, patch_size,
                     overlap, num_orientations, U, output, epsilon );

    return EXIT_SUCCESS;
}


void PLE_sub_routine(			
        MatrixXd * noisy
        ,			MatrixXd * invVar
        ,			int num_channels
        ,			int patch_size
        ,			int overlap
        ,			int num_orientations
        ,			MatrixXd * U
        ,			char const * output
        ,			double epsilon
        ){
    //initialize orientations and DCT
    int num_models = num_orientations + 1;

    //oriented bases
    NSym * covs;
    covs = initialize_covs( num_models, patch_size );

    //initial means are set to zero
    NColV * zero_mus = new NColV [num_models];
    int data_size = patch_size * patch_size;
    NColV mu( data_size );
    mu = 0;
    for( int i = 0; i < num_models; i++ )
        zero_mus[i] = mu;

    NColV * new_mus = new NColV[ num_models ];
    NSym * new_covs = new NSym[ num_models ];
    for( int i = 0; i < num_models; i++ ){
        new_mus[i] = zero_mus[i];
        new_covs[i] = covs[i];
    }

    MatrixXd * restored = PLE_denoise( noisy, invVar, overlap,
                                       patch_size, num_channels,
                                       num_models, U, new_covs, new_mus,
                                       epsilon);

    cout << "Writting final image." << endl;
    int num_channels_out = 1;
    MatrixXd * outToWrite = remake_bayer( restored, num_channels, num_channels_out );
    imwrite_exr( output, outToWrite, num_channels_out );

    delete [] zero_mus;
    delete [] covs;
    delete [] noisy;
    delete [] restored;
}

//no desire to optimize PLE because of its poor performance 
//hence an ugly interface strides between the two libraries
void convertFormat( 
        VectorXd **& eigen
        ,			NMatrix *& newmat
        ,			int map_rows
        ,			int map_cols
        ,			int data_size
        ,			bool inverse
        ){
    if( !inverse ){
        //VectorXd ** to NMatrix *
        newmat = new NMatrix [map_rows];
        //allocate and convert
        for(int row = 0; row < map_rows; row++){
            MatrixXd vmatrix( data_size, map_cols );
            vectorArray2Matrix( eigen, row, vmatrix );
            NMatrix temp( map_cols, data_size );
            temp << vmatrix.data();
            newmat[row] = temp.t();
            delete [] eigen[row];
        }
        //release
        delete[] eigen;
    }else{
        //NMatrix * to VectorXd **
        eigen = new VectorXd * [map_rows];
        //allocate and convert
        for(int row = 0; row < map_rows; row++){
            eigen[row] = new VectorXd [map_cols];
            for(int col = 0; col < map_cols; col++){
                VectorXd temp(data_size);
                NColV another = newmat[row].Column(col+1);
                for(int i = 0; i < data_size; i++)
                    temp(i) = another(i+1);
                eigen[row][col] = temp;
            }
        }
        //release
        delete [] newmat;
    }
}

//cluster + restore + update
MatrixXd * PLE_denoise( 	
        MatrixXd * noisy
        ,		MatrixXd * invVar
        ,		int overlap
        ,		int patch_size
        ,		int num_channels
        ,		int num_models
        ,		MatrixXd * U
        ,		NSym * covs
        ,		NColV * mus
        ,		double epsilon
        ){
    cout << "PLE denoising starts..." << endl;

    //general parameters
    int image_rows = noisy[0].rows();
    int image_cols = noisy[0].cols();
    int map_rows = num_patches( image_rows, patch_size, overlap );
    int map_cols = num_patches( image_cols, patch_size, overlap );
    int data_size = pow(patch_size, 2)*num_channels;

    //convert (noisy_patches will be released in convertFormat)
    VectorXd ** noisy_patches = image2patches( noisy, image_rows, image_cols, num_channels, overlap, patch_size );
    NMatrix * patches_at_row = NULL;
    convertFormat( noisy_patches, patches_at_row, map_rows, map_cols, data_size, false );
    for (int ii=0; ii<map_rows; ii++)
        for (int kk=1; kk<=data_size; kk++)
            for (int jj=1; jj<=map_cols; jj++)
                if (!isfinite(patches_at_row[ii](kk,jj))) {
                    cout << "Warning nan value in patches_at_row!!!! " << endl;
                    exit (EXIT_FAILURE);
                }


    //convert (U matrix transformation will be released in convertFormat)
    VectorXd ** U_patches = image2patches( U, image_rows, image_cols, num_channels, overlap, patch_size );
    NMatrix * U_at_row = NULL;
    convertFormat( U_patches, U_at_row, map_rows, map_cols, data_size, false );

    //convert (invVar matrix transformation will be released in convertFormat)
    VectorXd ** invVar_patches = image2patches( invVar, image_rows, image_cols, num_channels, overlap, patch_size );
    NMatrix * invVar_at_row = NULL;
    convertFormat( invVar_patches, invVar_at_row, map_rows, map_cols, data_size, false );

    //the main cycle
    MatrixXi cluster_map( map_rows, map_cols );

    //First iteration is treated separately to deal with color
    NMatrix * restored_patches = NULL;
    make_clusters_init( mus, covs, patches_at_row, invVar_at_row,
                        cluster_map, num_models, num_channels,
                        U_at_row, restored_patches );

    VectorXd ** restored = NULL;
    convertFormat( restored, restored_patches, map_rows, map_cols, data_size, true );
    MatrixXd * filtered = patches2image( restored, overlap, patch_size, num_channels, image_rows, image_cols, true );

    for( int row = 0; row < map_rows; row++ )
        delete [] restored[row];
    delete [] restored;

    //release and return
    delete [] patches_at_row;
    delete [] mus;
    delete [] covs;

    return filtered;
}


//classify a bunch of patches into K classes
//and each patch's label will be stored in cluster_map
void make_clusters_init( 	
        NColV const * mus
        ,		NSym const * covs
        ,		NMatrix const * noisy_patches
        ,		NMatrix * invVar_at_row
        ,		MatrixXi & cluster_map
        ,		int num_models
        ,		int num_channels
        ,		NMatrix * U
        ,		NMatrix *& restored_patches
        ){

    //	general parameters
    int data_size = noisy_patches[ 0 ].Nrows();
    int data_sizePatch = noisy_patches[ 0 ].Nrows() / num_channels;
    int map_rows = cluster_map.rows();
    int map_cols = cluster_map.cols();

    //	allocate memory for restored_patches
    if( restored_patches == NULL ){
        NMatrix temp( data_size, map_cols );
        restored_patches = new NMatrix [ map_rows ];
        for( int row = 0; row < map_rows; row++ )
            restored_patches[ row ] = temp;
    }

    //	inverse here rather than in the Wiener filter
    NSym * inv_cov = new NSym[num_models];
    double * cov_determinant = new double [num_models];

    for( int model = 0; model < num_models; model++ ){
        NMatrix tmp = covs[model].i();

        //inverse covariance
        NSym container(data_sizePatch);
        container << tmp;
        inv_cov[model] = container;

        //determinant
        cov_determinant[model] = log(covs[model].Determinant());
    }


    //	restore and cluster
#pragma omp parallel for schedule( static )
    for( int row  = 0; row < map_rows; row++ )
        for( int col = 1; col <= map_cols; col++ ){

            NColV Uv = U[ row ].Column( col );
            NColV invVarV = invVar_at_row[ row ].Column( col );
            NColV obsV = noisy_patches[ row ].Column( col );
            for ( int i=1; i<=data_sizePatch; i++ )
                if (!isfinite(obsV(i))) {
                    cout << "Warning nan value in make_clusters_init: OBSV !!!! " << endl;
                    exit (EXIT_FAILURE);
                }
            double likelihoods [ num_models ];
            NColV restored_aux(data_size);

            //	inverse here rather than in the Wiener filter
            NMatrix * Wiener_shrinkage = new NMatrix [num_models];


            for ( int ii=0; ii < num_models; ii++)
                likelihoods[ii] = 0.0;

            for ( int ch=0; ch < num_channels; ch++ ) {

                // get matrix U corresponding to current patch position
                // and transform it into a diagonal matrix
                NDiag u(data_sizePatch);
                for ( int i=1; i<=data_sizePatch; i++ )
                    u(i) = Uv(ch*data_sizePatch+i);

                // get matrix invVar corresponding to current patch position
                // and transform it into a diagonal matrix
                NDiag iVar(data_sizePatch);
                for ( int i=1; i<=data_sizePatch; i++ )
                    iVar(i) = invVarV(ch*data_sizePatch+i)*u(i);

                NDiag UtinvVarU = u.t()*iVar*u;
                NDiag UtinvVar = u.t()*iVar;

                //for( int model = 0; model < num_models; model++ )
                //Wiener_shrinkage[model] = (( UtinvVarU + inv_cov[model]).i())*UtinvVar;

                //NColV best_inpaint( data_size );
                NColV obs(data_sizePatch);
                for ( int i=1; i<=data_sizePatch; i++ )
                    obs(i) = obsV(ch*data_sizePatch+i);

                //				NColV obs = noisy_patches[ row ].Column( col ).segment(ch*data_sizePatch,data_sizePatch);

                //#pragma omp parallel for schedule(static)
                for( int model = num_models - 1; model >= 0; model-- )	{

                    Wiener_shrinkage[model] = (( UtinvVarU + inv_cov[model]).i())*UtinvVar;

                    //Wiener filtering
                    NColV observation_centered = obs - u*mus[model];
                    //best_guess = Wiener_shrinkage[model]*observation_centered;

                    //calculate the criterion
                    NColV estimate_centered = Wiener_shrinkage[model]*observation_centered;
                    NColV tmp = u*estimate_centered - observation_centered;
                    double likelihood = DotProduct( tmp, iVar * tmp );
                    likelihood += cov_determinant[model];
                    likelihood += DotProduct( estimate_centered, inv_cov[model] * estimate_centered );

                    likelihoods[ model ] += -1 * likelihood;
                }
            } // end for channels

            //if there is no clear winner, let's assign the patch to DCT
            double highest_likelihood = likelihoods[ num_models - 1 ];
            int winner = num_models - 1;
            for( int model = num_models - 2; model >= 0; model-- )
                if( likelihoods[ model ] > highest_likelihood ){
                    highest_likelihood = likelihoods[ model ];
                    winner = model;
                }

            for ( int ch=0; ch < num_channels; ch++ ) {

                // get matrix U corresponding to current patch position
                // and transform it into a diagonal matrix
                NDiag u(data_sizePatch);
                for ( int i=1; i<=data_sizePatch; i++ )
                    u(i) = Uv(ch*data_sizePatch+i);

                // get matrix invVar corresponding to current patch position
                // and transform it into a diagonal matrix
                NDiag iVar(data_sizePatch);
                for ( int i=1; i<=data_sizePatch; i++ )
                    iVar(i) = u(i)*invVarV(ch*data_sizePatch+i);

                NDiag UtinvVarU = u.t()*iVar*u;
                NDiag UtinvVar = u.t()*iVar;

                NMatrix	Wiener_shrinkageChosen = (( UtinvVarU + inv_cov[winner]).i())*UtinvVar;

                NColV obs(data_sizePatch);
                for ( int i=1; i<=data_sizePatch; i++ )
                    obs(i) = obsV(ch*data_sizePatch+i);

                //Wiener filtering
                NColV observation_centered = obs - u*mus[winner];
                NColV best_guess = Wiener_shrinkageChosen*observation_centered;

                for ( int i=1; i<=data_sizePatch; i++ )
                    if (!isfinite(obs(i))) {
                        cout << "Warning nan value in make_clusters_init: OBS !!!! " << endl;
                        exit (EXIT_FAILURE);
                    }

                for ( int i=1; i<=data_sizePatch; i++ )
                    if (isfinite(best_guess(i)))
                        restored_aux(ch*data_sizePatch + i) = best_guess(i);
                    else {
                        cout << "Warning nan value in make_clusters_init !!!! " << endl;
                        exit (EXIT_FAILURE);
                    }
            }

            restored_patches[ row ].Column(col) = restored_aux;
            cluster_map( row, col - 1 ) = winner;

            //	release
            delete [] Wiener_shrinkage;

        } // end for patches

    delete [] inv_cov;
    delete [] cov_determinant;
    //delete [] Wiener_shrinkage;

}

int most_probable_generator( 	
        NColV const * mus
        ,		double const * cov_determinant
        ,		NSym const * inv_cov
        ,		NMatrix const * Wiener_shrinkage
        ,		double sigma
        ,		int num_models
        ,		NDiag U
        ,		NDiag iVar
        ,		NColV const & observation
        ,		NColV & guessed_signal
        ,		int row
        ,		int col
        ){

    //guessed_array contains the estimates from various model aassumptions
    NColV * guessed_array = new NColV [ num_models ];
    double likelihoods [ num_models ];

#pragma omp parallel for schedule(static)
    for( int model = num_models - 1; model >= 0; model-- )
        likelihoods[ model ] = guess_signal( mus[model], cov_determinant[model], inv_cov[model], Wiener_shrinkage[model], sigma, U, iVar, observation, guessed_array[model], row, col, model );


    //if there is no clear winner, let's assign the patch to DCT
    double highest_likelihood = likelihoods[ num_models - 1 ];
    int winner = num_models - 1;
    for( int model = num_models - 2; model >= 0; model-- )
        if( likelihoods[ model ] > highest_likelihood ){
            highest_likelihood = likelihoods[ model ];
            winner = model;
        }

    //final verdict
    guessed_signal = guessed_array[ winner ];
    delete [] guessed_array;
    return winner;
}

//Wiener filtering, not optimized
double guess_signal( 	
        NColV const & mu
        ,		double cov_determinant
        ,		NSym const & inv_cov
        ,		NMatrix const & shrinkage
        ,		double sigma
        ,		NDiag U
        ,		NDiag iVar
        ,		NColV const & observation
        ,		NColV & best_guess
        ,		int row
        ,		int col
        ,		int model
        ){
    int data_size = mu.Nrows();
    double squared_sigma = pow( sigma, 2 );

    //Wiener filtering
    NColV observation_centered = observation - U*mu;
    best_guess = shrinkage*observation_centered;

    //calculate the criterion
    NColV estimate_centered =  best_guess;
    NColV tmp = U*estimate_centered - observation_centered;
    double likelihood = DotProduct( tmp, iVar * tmp );
    likelihood += cov_determinant;
    likelihood += DotProduct( estimate_centered, inv_cov * estimate_centered );

    best_guess = best_guess + mu;

    //use the highest in the next step, so add a minus sign here
    return -1 * likelihood;
}

//initialize Gaussian models
NSym * initialize_covs( 
        int num_models
        ,		int patch_size
        ){

    //	allocate memory
    int data_size = pow( patch_size, 2 );
    NSym temp( data_size );
    NSym * covs = new NSym [ num_models ];
    for( int model = 0; model < num_models; model++ )
        covs[ model ] = temp;

    //	set eigenvalues
    NDiag D( data_size );
    bool use_HS = false;  //heavyside
    for( int i = 1; i <= data_size; i++ )
        if( !use_HS )
            //			adjust the decay speed, or it won't work for large patch_size: really a hack
            D(i) = pow(2, 20.5-0.5/(exp((patch_size-8)/5))*(double)i) +  0.1;
        else{
            if( i < 6 )
                D(i) = 1e4;
            else
                D(i) = 1e-1;
        }

    //	set the eigenvectors for the oriented eigenvectors
    int num_orientations = num_models - 1;
#pragma omp parallel for
    for( int model = 0; model < num_orientations; model++ ){
        NMatrix U( data_size, data_size );
        int n_samples = 3 * data_size * data_size;
        PCA( M_PI*(double)model/(double)num_orientations, n_samples, patch_size, U );
        NMatrix cov = U*D*U.t();
        covs[ model ] << cov;
    }

    //	DCT
    NMatrix U( data_size, data_size );
    get_DCT_basis(U);
    covs[ num_orientations ] << U*D*U.t();

    return covs;
}



//create DCT to complete the oriented bases
void get_DCT_basis( NMatrix & output ){
    //	the output matrix must be of dimension ( n*n, n*n )
    //	n is the matrix's dimension
    int n = (int) sqrt(output.Ncols());
    NMatrix basis(n,n);
    //	first build the one dimensional basis
    int k = (n - n%2)/2;
    int col_index = 1;
    NColV container(n);
    //	the component frequency in arranged in ascending order
    for( int i = 0; i <= k; i++ ){
        //		cos is always good
        for( int counter = 0; counter < n; counter++ )
            container(counter+1) = cos((2.*(double)M_PI*i/(double)n*counter));
        container /= sqrt(container.SumSquare());
        basis.Column(col_index) << container;
        col_index++;
        //		sin needs some care
        for( int counter = 0; counter < n; counter++ )
            container(counter+1) = sin((2.*(double)M_PI*i/(double)n*counter));
        double norm_square = container.SumSquare();
        if( norm_square > 1e-4 ){
            basis.Column(col_index) << container/sqrt(norm_square);
            col_index++;
        }
    }
    NMatrix product(n,n);
    NColV another(n);
    col_index = 1;
    //	try to order the basis vectors so that their frequency only increase
    for( int counter = 2; counter <= 2*n; counter ++ )
        for( int i = 1; i < n+1; i++ ){
            container << basis.Column(i);
            int j = counter - i;
            if( j > n || j < 1 )
                continue;
            another << basis.Column(j);
            product << container*another.t();
            //			reorder the vectors according to (i+j)
            output.Column(col_index) = product.AsMatrix(n*n,1);
            col_index += 1;
        }
    cout << "got the DCT basis, a check here (result shall be " << n*n << " ): " <<  (output*output.t()).Sum() << endl;
}

//perform PCA for n_samples patches drawn from synthetic images 
void PCA( 	
        double theta
        ,		int n_samples
        , 		int patch_size
        ,		NMatrix & U
        ){

    int data_size = pow( patch_size, 2 );
    NMatrix result( data_size, n_samples );
    //sample
    get_synthesized_patches( theta, n_samples, patch_size, result );

    //calc the mean
    NColV mean( data_size );
    mean = 0;
    for( int sample = 1; sample <= n_samples; sample++ )
        mean += result.Column( sample );
    mean /= (double) n_samples;

    //calc the cov
    NMatrix C( data_size, data_size );
    C = 0;
    for( int sample = 1; sample <= n_samples; sample++ ){
        NColV container = result.Column( sample ) - mean;
        C += container * container.t();
    }
    C /= (double) n_samples;
    NSym cov( data_size );
    cov << C;

    //PCA by SVD
    NDiag D( data_size );
    SVD( cov, D, U );
}

//get n_samples synthetic patches from some feature images 
void get_synthesized_patches( 	
        double theta
        ,				int n_samples
        ,		 		int patch_size
        ,		 		NMatrix & result
        ){
    static int starter = 0;
    if( starter == 0 ){
        init_randmt_auto();
        starter = 1;
    }

    cout << "synthesize the image feature orientation = " << theta * 180 / (double) M_PI << endl;

    //get the image not yet blurred by gaussian kernels
    int image_dim = 100, kernel_dim = 30;
    NMatrix image( image_dim, image_dim );
    get_synthesized_image( theta, image );

    //blur it
    NMatrix kernel( kernel_dim, kernel_dim );
    int valid_dim = image_dim - kernel_dim + 1;
    NMatrix blurred_image( valid_dim, valid_dim );
    //just an arbitrary number
    int num_kernels = 4;
    NMatrix * images = new NMatrix [ num_kernels ];
    double sigma = 0;
    for( int i = 1; i <= num_kernels; i++ ){
        sigma = (double)i * 2;
        Gaussian_kernel( sigma, kernel );
        convolute( image, kernel, blurred_image );
        images[ i - 1 ] = blurred_image;
    }


    //sampling
    int data_size = pow( patch_size, 2 );
    //image_index
    int index;
    //some parameters to improve sampling efficiency
    double r, factor;
    if( theta <= (double) M_PI/4 || theta >= (double) M_PI/4*3 )
        factor = max( 1/sqrt((double)2.), abs(cos(theta)) );
    else
        factor = max( 1/sqrt((double)2.), abs(sin(theta)));
    //r is the maximal distance a sample patch's can go w.r.t. the image's center.
    r = (double)(valid_dim-patch_size)/(factor*2);
    double center_x, center_y;
    double rho, whichSide;
    int ncenter_x, ncenter_y, left_j, right_j, upper_i, lower_i;
    for( int i = 1; i <= n_samples; i++ ){
        //first decide the image to sample from
        index = (int)(rand_unif() * (double)num_kernels);
        blurred_image = images[index];
        while(true){
            //image's center
            center_x = (1+valid_dim)/2;
            center_y = center_x;
            //determine the sample coordinates
            rho = r*rand_unif();
            whichSide = ((int)(rand_unif() + 0.5) - 0.5)*2;
            //sample along the gradient
            ncenter_y = (int)(center_y + rho*whichSide*cos(theta));
            ncenter_x = (int)(center_x - rho*whichSide*sin(theta));
            //more precisely
            left_j = ncenter_x - (int)(patch_size/2);
            right_j = left_j + patch_size - 1;
            upper_i = ncenter_y - (int)(patch_size/2);
            lower_i = upper_i + patch_size - 1;
            //check if the patch lies in the image
            if( (left_j >= 1) && (right_j <= valid_dim) && (upper_i >= 1) && (lower_i <= valid_dim))
                break;
        }
        //scrap useless (i.e. uniform) pattern types since they have nothing to do with the feature
        NMatrix patch = blurred_image.SubMatrix( upper_i, lower_i, left_j, right_j );
        double var = patch.SumSquare()/(double)data_size - pow( patch.Sum()/(double)( data_size ), 2 );
        if( var < 1e-4 ){
            i = i - 1;
            continue;
        }
        //format the patch for storage
        result.Column( i ) << patch.AsMatrix( data_size, 1 );
    }
    delete [] images;
}


//produce a synthetic image with specified feature gradient 
void get_synthesized_image( 	
        double theta
        ,			NMatrix & image
        ){
    int n = image.Nrows();
    //make sure that the black/white borderline goes through the image's center
    double center_x = (double) (n+1)/2;
    double center_y = (double) (n+1)/2;
    for(double row = 1; row <= n; row++)
        for(double col = 1; col <= n; col++){
            //cope with the singularity
            if( abs(theta - (double) M_PI/2) < 1e-4 ){
                if( col <= center_x )
                    image( row, col ) = 1;//0;
                else
                    //image( row, col ) = 1;
                    image( row, col ) = 0.5;
            }
            else{
                //normal scenario
                if( tan(theta)*(col - center_x) > row - center_y )
                    //image( row, col ) = 1;
                    image( row, col ) = 1;
                else
                    image( row, col ) = 0.5;//0;
            }
        }
}

//built on a n*n matrix and centered on (n+1/2,n+1/2)
//with its standard deviation specified by sigma
//the whole thing is then stored to the NMatrix called kernel
void Gaussian_kernel( 
        double sigma
        ,			NMatrix & kernel
        ){
    int n = kernel.Nrows();
    //mean
    NColV mu(2);
    mu = 0;
    //cov
    NMatrix cov( 2, 2 );
    cov = 0;
    for( int row = 1; row <= 2; row++ )
        cov( row , row ) = pow( sigma, 2 );
    NSym sym_cov( 2 );
    sym_cov << cov;
    //fill in the kernel
    for( int row = 1; row <= n; row++ )
        for(int col = 1; col <= n; col++ ){
            NColV x( 2 );
            x << row  << col;
            //centering
            x -= (double)(n+1)/2.;
            kernel( row , col ) = gaussian_density( x, mu, sym_cov );
        }
    //normalize to make it a discrete probability
    double total = kernel.Sum();
    kernel /= total;
}


//self evident
double gaussian_density(	
        NColV const & x
        ,		NColV const & mu
        ,		NSym const & cov
        ){
    double data_size = (double) mu.Nrows();
    NEWMAT::LowerTriangularMatrix L = Cholesky( cov );
    double determinant = 1;
    for( int i = 1; i <= data_size; i++ )
        determinant *= pow( L( i,i ), 2 );
    double denominator = pow( determinant, 0.5 );
    NColV y = L.i()*(x-mu);
    double numerator = exp(-0.5*y.SumSquare());
    return numerator/denominator;
}

//convolute a synthesized feature image depicting
//a purely black region borders on a purely white one
void convolute( 
        NMatrix const & image
        ,		NMatrix const & kernel
        ,		NMatrix & result
        ){
    int result_dim = result.Nrows();
    int kernel_dim = kernel.Nrows();
    for( int row = 1; row <= result_dim; row++ )
        for( int col = 1; col <= result_dim; col++ ){
            NMatrix mask = image;
            mask = 0;
            mask.SubMatrix( row, row+kernel_dim-1, col, col+kernel_dim-1 ) = kernel;
            //elementwise Schur Product
            result(row,col) = SP(mask, image).Sum();
        }
}

void write_covs( const char *fname, NSym * covs, int num_models, int patch_size ) {

    int data_sz = patch_size*patch_size;
    //FILE *fp = fopen("covs.txt","w");
    FILE *fp = fopen(fname,"w");

    for (int c=0; c < num_models; c++) {

        for (int i=1; i<= data_sz; i++)
            for (int j=1; j<= data_sz; j++)
                fprintf(fp,"%f ", covs[c](i,j));
    }

    fclose(fp);
}

void write_mus( const char *fname, NColV * mus, int num_models, int patch_size ) {

    int data_sz = patch_size*patch_size;
    //FILE *fp = fopen("mus.txt","w");
    FILE *fp = fopen(fname,"w");

    for (int c=0; c < num_models; c++) {

        for (int i=1; i<= data_sz; i++)
            fprintf(fp,"%f ", mus[c](i));
    }

    fclose(fp);
}


NSym *read_covs( const char *fname, int num_models, int patch_size ) {

    int data_size = patch_size*patch_size;
    NSym temp( data_size );
    NSym * covs = new NSym [ num_models ];
    for( int model = 0; model < num_models; model++ )
        covs[ model ] = temp;

    FILE *fp = fopen(fname,"r");
    if ( !fp )
        cout << "File covs.txt coud not been read." << endl;
    //double *data = (double*) malloc(sizeof(double)*num_models*data_size*data_size);
    double *data = (double*) malloc(sizeof(double));

    for (int c=0; c < num_models; c++)
        for (int i=1; i<= data_size; i++)
            for (int j=1; j<= data_size; j++) {
                int readVal = fscanf(fp,"%lf ",data);
                if (readVal == 1)
                    covs[c](i,j) = *data;
                else {
                    cout << readVal << " values written from the covs.txt file." << endl;
                    exit (EXIT_FAILURE);
                }
            }

    fclose(fp);

    return covs;
}

NColV *read_mus( const char *fname, int num_models, int patch_size ) {

    int data_size = patch_size * patch_size;
    NColV * mus = new NColV [num_models];
    NColV mu( data_size );
    mu = 0;
    for( int i = 0; i < num_models; i++ )
        mus[i] = mu;

    FILE *fp = fopen(fname,"r");
    if ( !fp )
        cout << "File mus.txt coud not been read." << endl;
    double *data = (double*) malloc(sizeof(double));
    for (int c=0; c < num_models; c++)
        for (int i=1; i<= data_size; i++)	{
            int readVal = fscanf(fp,"%lf ",data);
            if ( readVal == 1 )
                mus[c](i) = *data;
            else {
                cout << readVal << " values written from the mus.txt file." << endl;
                exit (EXIT_FAILURE);
            }
        }

    fclose(fp);

    return mus;
}








