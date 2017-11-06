/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file hbe.cpp
 * @brief HBE denoising functions
 *
 * @author Cecilia Aguerrebere adapted from Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <float.h>

#include "hbe.h"
#include "LibMatrix.h"
#include "LibImages.h"
#include "Utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

void initializeParameters(
	nlbParams &o_paramStep2
,   const double p_sigma
,	const ImageSize &p_imSize
,   const bool p_verbose
,	unsigned pSize
,	unsigned offset
){
	//! Standard deviation of the noise
	o_paramStep2.sigma = p_sigma;

	//! Size of patches
	if (p_imSize.nChannels == 1) {				
		o_paramStep2.sizePatch = pSize;		
	}
	else {
        o_paramStep2.sizePatch = pSize;
	}

	//! Number of similar patches
	if (p_imSize.nChannels == 1) {
		o_paramStep2.nSimilarPatches =	(p_sigma < 20. ? 2*pSize*pSize :
										(p_sigma < 40. ? 25 :
										(p_sigma < 80. ? 30 : 45)));
	}
	else {
		o_paramStep2.nSimilarPatches = o_paramStep2.sizePatch * o_paramStep2.sizePatch * 3;
	}

	//! Offset: step between two similar patches
	o_paramStep2.offSet = offset;

	//! Size of the search window around the reference patch (must be odd)	
	o_paramStep2.sizeSearchWindow = pSize*3+1;
	if (o_paramStep2.sizeSearchWindow % 2 == 0) {
		o_paramStep2.sizeSearchWindow++;
	}

	//! Size of boundaries used during the sub division
	o_paramStep2.boundary = int(1.5 * double(o_paramStep2.sizeSearchWindow));

	//! Print information?
	o_paramStep2.verbose = p_verbose;

}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy image;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_imSize: size of the image;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/
int runHBE(
	std::vector<double> const& imNoisy
,   std::vector<double> &imBasic
,	std::vector<double> &imFinal
,   std::vector<double> &imUmask
,	const ImageSize &p_imSize
,	const double p_sigma
,   const bool p_verbose
,	double alphaH
,	double alphaL
,	const double minPixKnown
,	unsigned pSize
,	unsigned offset
,	double Nfactor
,	double NfactorPrior
,	double epsilon_pd
,	double varSigma
){
	//! Only 1, 3 or 4-channels images can be processed.
	const unsigned chnls = p_imSize.nChannels;
	if (! (chnls == 1 || chnls == 3 || chnls == 4)) {
		cout << "Wrong number of channels. Must be 1 or 3!!" << endl;
		return EXIT_FAILURE;
	}

	//! Number of available cores
	unsigned nbThreads = 1;
#ifdef _OPENMP
    nbThreads = omp_get_max_threads();
    if (p_verbose) {
        cout << "Open MP is used" << endl;
    }
#endif
	
	//! Initialize output image
	imFinal.resize(imNoisy.size());
	for (unsigned i=0; i<p_imSize.whc; i++ )
		imFinal[i] = 0.0;

	//! Parameters Initialization
	nlbParams paramStep2;
	initializeParameters(paramStep2, p_sigma, p_imSize, 
							p_verbose, pSize, offset);

	//! Start processing
	if (paramStep2.verbose) {
		cout << "Start processing...";
		cout << "Patch size: " << paramStep2.sizePatch << endl;
		cout << "Num similar patches: "  << paramStep2.nSimilarPatches << endl;		
		cout << "Search window size: "  << paramStep2.sizeSearchWindow << endl;
		cout << "Boundary: "  << paramStep2.boundary << endl;
		cout << "Offset: "  << paramStep2.offSet << endl;
		cout << " " << endl;
	}

	//! Divide the noisy image into sub-images in order to easier 
	//!parallelize the process
	const unsigned nbParts = 2 * nbThreads;
	vector<vector<double> > imNoisySub(nbParts), imUmaskSub(nbParts), 
	imPLEmodelSub(nbParts), imBasicSub(nbParts), imFinalSub(nbParts);
	ImageSize imSizeSub;
	cout << "NbParts: " << nbParts << endl;	
	
	unsigned maxIts = 3;
	//! Main outer loop
	for (unsigned it=0; it<maxIts; it++ ) {
	
		if ( it > 0 ) 
			for (unsigned k = 0; k < p_imSize.whc; k++) {
				imBasic[k] = imFinal[k];
				imFinal[k] = 0.0;
			}	
	
		if (subDivide(imNoisy, imNoisySub, p_imSize, imSizeSub, 
						paramStep2.boundary, nbParts)
			!= EXIT_SUCCESS) {
			return EXIT_FAILURE;
		}
		if (subDivide(imBasic, imBasicSub, p_imSize, imSizeSub, 
						paramStep2.boundary, nbParts)
			!= EXIT_SUCCESS) {
			return EXIT_FAILURE;
		}
		if (subDivide(imUmask, imUmaskSub, p_imSize, imSizeSub, 
					paramStep2.boundary, nbParts)
			!= EXIT_SUCCESS) {
			return EXIT_FAILURE;
		}
		if (subDivide(imFinal, imFinalSub, p_imSize, imSizeSub, 
					paramStep2.boundary, nbParts)
			!= EXIT_SUCCESS) {
			return EXIT_FAILURE;
		}
		
		//! Process all sub-images
	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic, nbParts/nbThreads) \
				shared(imNoisySub, imBasicSub, imFinalSub) \
			firstprivate (paramStep2)
	#endif
		for (int n = 0; n < (int) nbParts; n++) {
			processHBE(imNoisySub[n], imBasicSub[n], imUmaskSub[n], 
							imFinalSub[n], imSizeSub, paramStep2, alphaH, 
							alphaL, minPixKnown, Nfactor, 
							NfactorPrior, epsilon_pd, varSigma);
		}

		//! Get the final result
		if (subBuild(imFinal, imFinalSub, p_imSize, imSizeSub, paramStep2.boundary)
			!= EXIT_SUCCESS) {
			return EXIT_FAILURE;
		}

		if (paramStep2.verbose) {
			cout << "done." << endl << endl;
		}
	}
	
	return EXIT_SUCCESS;
}

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_imSize: size of i_imNoisy;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
void processHBE(
	std::vector<double> const& imNoisy
,	std::vector<double> &imBasic
,	std::vector<double> &Umask
,	std::vector<double> &imFinal
,	const ImageSize &p_imSize
,	nlbParams &p_params
,	double alphaH
,	double alphaL
,   double minPixKnown
,	double Nfactor
,	double NfactorPrior
,	double epsilon_pd
,	double varSigma
){
	//! Parameters initialization
	const unsigned sW		= p_params.sizeSearchWindow;
	const unsigned sP		= p_params.sizePatch;
	const unsigned sP2		= sP * sP;
	const unsigned sPC		= sP2 * p_imSize.nChannels;
	const unsigned nSP		= p_params.nSimilarPatches;
	double minPixKnownP = minPixKnown*sPC;
	
	vector<double> Wout;
	cout << "Allocating memory for Wout" << endl;
	Wout.resize(sP2*sP2*300);	
	cout << "End allocating memory for Wout" << endl;

	//! Allocate Sizes
	imFinal.resize(p_imSize.whc);

	//! Used matrices during estimate
	vector<vector<double> > group3d(p_imSize.nChannels, vector<double> (nSP * sP2));
	vector<double> group3dNoisy(sW * sW * sPC), group3dBasic(sW * sW * sPC),
	group3dUmask(sW * sW * sPC),group3dBasicPrior(sW * sW * sPC);
	vector<unsigned> index(sW * sW);

	//! Mask: non-already processed patches
	vector<bool> mask(p_imSize.wh, false);

	//! Only pixels of the center of the image must be processed (not the boundaries)
	for (unsigned i = sW; i < p_imSize.height - sW; i++) {
		for (unsigned j = sW; j < p_imSize.width - sW; j++) {
			mask[i * p_imSize.width + j] = true;
		}
	}
	
	//! ponderation: weight sum per pixel
	vector<double> weight(imNoisy.size(), 0.);
	vector<double> Umask_current = Umask;
	
	double nu2_curr = alphaL;
	double kappa2_curr = alphaL;
	
	unsigned totalPixToFill = (p_imSize.height-2*sW) * (p_imSize.width-2*sW);	
	unsigned nSimP, nSimPprior;
	unsigned percentageToFill = 100;
	
	NDiag Ew(sPC);
	for (unsigned k=1; k<=sPC; k++)
		Ew(k) = varSigma;
	
	//! For all patches...	
	for (unsigned ij = 0; ij < p_imSize.wh; ij += p_params.offSet) {	
		
		//! Variables to keep track of restoration status
		unsigned old_percentage = percentageToFill;
		unsigned pixToFill = pixels_to_fillMask( mask, p_imSize.wh);
		percentageToFill = 100 * pixToFill / totalPixToFill;		
		if (( (int)percentageToFill % 20 == 0) && (percentageToFill < old_percentage) )
			cout << "Still " << percentageToFill << 
			"% of the pixels to be filled (" << pixToFill <<
			" of " << totalPixToFill << ")."<< endl;
		
		//! Only non-seen patches are processed
		if ( mask[ij] ) {

			//! Search for similar patches around the reference one
			nSimP = p_params.nSimilarPatches;
			nSimP = findSimilarPatches(imNoisy, imBasic, Umask, 
											group3dNoisy, group3dBasic, 
											group3dBasicPrior, group3dUmask, 
											index, ij, p_imSize, p_params, 
											Nfactor, NfactorPrior, 
											minPixKnownP, &nSimPprior);	

			//! Count number of known pixels in the current patch
			double knownPix = 0.0;
			for (unsigned p = 0; p < sP; p++) 
				for (unsigned q = 0; q < sP; q++) 
					knownPix += Umask[ij + p * p_imSize.width + q];												

			//! Depending on number of known pixels and similar patches,
			//! set the current values of kappa and nu
			if ( (nSimP < 10) || (knownPix < minPixKnownP))	{							
				nu2_curr = alphaH;
				kappa2_curr = alphaH;								
			} else {
				nu2_curr = alphaL;
				kappa2_curr = alphaL;
			}

			//! Compute the restored patch
			//cout << "Restoring " << nSimP << " sim patch" << endl;
			computeBayesEstimate(group3dNoisy, group3dBasic, group3dBasicPrior, 
								group3dUmask, p_imSize, p_params, nSimP, 
								nSimPprior, kappa2_curr, nu2_curr,
								epsilon_pd, Ew, Wout);
			//cout << "End restoring " << nSimP << " sim patch" << endl;
			//! Perform aggregation
			computeAggregation(imFinal, weight, mask, group3dBasic, index, 
								p_imSize, p_params, nSimP);										
		}
		
	} //! End for ij

	//! Weighted aggregation
	computeWeightedAggregation(imNoisy, imBasic, imFinal, weight, p_params, p_imSize);

}


unsigned findSimilarPatches(
	std::vector<double> const& imNoisy
,	std::vector<double> const& imBasic
,	std::vector<double> const& Umask
,	std::vector<double> &o_group3dNoisy
,	std::vector<double> &o_group3dBasic
,	std::vector<double> &o_group3dBasicPrior
,	std::vector<double> &o_Umask
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const ImageSize &p_imSize
,	const nlbParams &p_params
,   double N
,   double Nprior
,	double minPixKnown
,	unsigned *nSimPprior_out
){
	//! Initialization
	const unsigned width	= p_imSize.width;
	const unsigned chnls	= p_imSize.nChannels;
	const unsigned wh		= width * p_imSize.height;
	const unsigned sP		= p_params.sizePatch;
	const unsigned sW		= p_params.sizeSearchWindow;
	const unsigned ind		= p_ij - (sW - 1) * (width + 1) / 2;
	vector<pair<double, unsigned> > distance(sW * sW);	
	
	vector<unsigned> prior_index(sW * sW);
	
		//! Compute distance between patches in the search window sW x sW
		for (unsigned i = 0; i < sW; i++) {
			for (unsigned j = 0; j < sW; j++) {
				const unsigned k = i * width + j + ind;
				double diff = 0.0;								
				double weightTotal = 0.0;		
				//! Compute weighted distance between the current patch 
				//! (indp) and each patch of the search window (indk)	
				for (unsigned c = 0; c < chnls; c++) {
					const unsigned dc = c * wh;		
					unsigned pij_aux = dc + p_ij;
					unsigned k_aux = dc + k;		
					for (unsigned p = 0; p < sP; p++) {
						unsigned pij_aux1 = pij_aux + p * width;
						unsigned k_aux1 = k_aux + p * width;
						for (unsigned q = 0; q < sP; q++) {	
							unsigned indp = pij_aux1 + q;
							unsigned indk = k_aux1 + q;
							double knownPix = Umask[indp]*Umask[indk];						
							double weight = knownPix > 0.0 ? 1.0 : 0.01;																			
							const double tmpValue = imBasic[indp] - imBasic[indk];						
							diff += tmpValue * tmpValue * weight;
							weightTotal += weight;							
						}
					}
				}
				
				diff = diff / weightTotal;								
				//! Save all distances
				distance[i * sW + j] = make_pair(diff, k);
			}
		}
				
	//! Keep only the nSimilarPatches best similar patches
	//! Esto a la salida da ordenados en ascending order los elementos 
	//! primeros p_params.nSimilarPatches
	//! elementos. Los elementos que siguen estan desordenados.
	partial_sort(distance.begin(), distance.begin() + p_params.nSimilarPatches, 
				distance.end(), comparaisonFirst);

	//! Define thresholds to select similar patches
	double threshold = N*distance[1].first;	
	if (threshold < DBL_MAX ) {
		unsigned n = 2;
		while ( (threshold <= 0) && (n < p_params.nSimilarPatches) 
					&& (distance[n].first < DBL_MAX)) {
			cout << "Threshold negativo " << threshold << endl;
			threshold = N*distance[n].first;
			n++;
		}	
	} else 
		threshold = 0.0;

	double threshold2 = Nprior*threshold/N > 
						distance[p_params.nSimilarPatches - 1].first 
						? Nprior*threshold/N : 
						distance[p_params.nSimilarPatches - 1].first;
						
	unsigned nSimP = 0;
	unsigned nSimPprior = 0;

	//! Register position of similar patches
	for (unsigned n = 0; n < distance.size(); n++) {
		if (distance[n].first <= threshold) {
			o_index[nSimP++] = distance[n].second;
		}
		if (distance[n].first <= threshold2) {
			prior_index[nSimPprior++] = distance[n].second;
		}
	}						

	//! Save similar patches into 3D groups
	for (unsigned c = 0, k = 0; c < chnls; c++) {
		unsigned cwh = c * wh;
		for (unsigned p = 0; p < sP; p++) {
			unsigned cwh_pwidth = cwh + p * width;
			for (unsigned q = 0; q < sP; q++) {
				unsigned cwh_pwidth_q = cwh_pwidth + q;
				
				for (unsigned n = 0; n < nSimP; n++, k++) {
					unsigned indk = o_index[n] + cwh_pwidth_q;
					o_group3dNoisy[k] = imNoisy[indk];
					o_group3dBasic[k] = imBasic[indk];					
					o_Umask[k] = Umask[indk];					
				}
	
			}
		}
	}
	

	//! Idem with the patches to compute prior
	for (unsigned c = 0, k = 0; c < chnls; c++) {
		unsigned cwh = c * wh;
		for (unsigned p = 0; p < sP; p++) {
			unsigned cwh_pwidth = cwh + p * width;
			for (unsigned q = 0; q < sP; q++) {
				unsigned cwh_pwidth_q = cwh_pwidth + q;
								
				for (unsigned n = 0; n < nSimPprior; n++, k++) {
					unsigned indk = prior_index[n] + cwh_pwidth_q;
					o_group3dBasicPrior[k] = imBasic[indk];
				}
			}
		}
	}

	*nSimPprior_out = nSimPprior;
	
	return nSimP;
}


/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_group3d: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *		- group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/


void computeBayesEstimate(
	std::vector<double> &i_group3dNoisy
,	std::vector<double> &io_group3dBasic
,	std::vector<double> &io_group3dBasicPrior
,	std::vector<double> &i_Umask
,	const ImageSize &p_imSize
,	nlbParams p_params
,	const unsigned p_nSimP
,	const unsigned p_nSimPprior
,	double kappa
,	double nu
,	double epsilon_pd
,	NDiag Ew
,	std::vector<double> &Wout
){
	//! Parameters initialization
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch * p_imSize.nChannels;

	//! Keep the non centered version of the data
	//! Patches are centered to comupte covariance matrices
	//! but the non-centered versions are also kept to compute
	//! the restoration
	vector<double> io_group3dBasicOrig;
	io_group3dBasicOrig = io_group3dBasic;
	vector<double> i_group3dNoisyOrig;
	i_group3dNoisyOrig = i_group3dNoisy;
	vector<double> i_UmaskBasic;
	i_UmaskBasic = i_Umask;

	//! Compute mu estimate for prior and CENTER data to compute the
	//! prior covariance matrix
	NColV mu_est_prior(sPC);
	mu_est_prior = centerDataForPrior(io_group3dBasicPrior, p_nSimPprior, sPC);
	
	//! Compute prior covariance matrix for prior
	NSym cov_est_prior = covarianceMatrixForPrior(io_group3dBasicPrior, 
							p_nSimPprior, sPC, epsilon_pd);	
	
	//! Make sure the zero pixels in the degradation mask are kept to zero in 
	//! the current patches
	for (unsigned i = 0; i < p_nSimP; i++) 				
		for (unsigned j = 0; j < sPC; j++) 						
			io_group3dBasic[j * p_nSimP + i] = i_Umask[j * p_nSimP + i]*
											   io_group3dBasic[j * p_nSimP + i];		


	//! Main loop that iterates between two steps:
	//! (i)  patch restoration and model mean computation
	//! (ii) model covariance matrix computation 
	NSym cov_est_basic = cov_est_prior;	
	unsigned maxIts = 2;
	
	for (unsigned iter = 0; iter < maxIts; iter++) {	
		//! Compute the filter W only once because it is the same 
		//! for all patches restored together
        computeW( i_group3dNoisyOrig, i_UmaskBasic, Wout,
                cov_est_basic, p_nSimP, sPC, p_nSimP, Ew );
		
		//! Compute the model mean	
		NColV mu_est_basic = centerData(i_group3dNoisyOrig, i_UmaskBasic, 
							 mu_est_prior, p_nSimP, sPC, kappa, p_nSimP, Wout);

		//! Restore the patches
		NColV f;
		NMatrix W(sPC,sPC);
		//! For each patch...
		for (unsigned p = 0; p < p_nSimP; p++) {
			NDiag U(sPC);
			U = 0.0;
			NColV noisy(sPC);			
			for (unsigned k = 1; k <= sPC; k++) {
				unsigned ind = (k-1)*p_nSimP + p;				
				noisy(k) = i_group3dNoisyOrig[ind];    //! value already centered
				U(k) = i_Umask[ind];      
			}

			unsigned k=0;
			for (unsigned ii = 1; ii <= sPC; ii++) 
				for (unsigned jj = 1; jj <= sPC; jj++,k++) 
					W(ii,jj) = Wout[sPC*sPC*p + k];				
					
			//! Compute restored patch (the mean will be added later)							
			f = W*noisy;
			
			//! Save restored CENTERED patch
			for (unsigned k = 0; k < sPC; k++)
				io_group3dBasic[k*p_nSimP + p] = f(k+1); 
		}

		//! Compute the covariance matrix of the set of similar patches		
		cov_est_basic = covarianceMatrix(io_group3dBasic, i_UmaskBasic, 
							mu_est_prior, cov_est_prior, 
							p_nSimP, sPC, kappa, nu, mu_est_basic, 
							p_params.beta, p_nSimP, 
							epsilon_pd, Ew);			

		//! Add mean to previously restored patches
		for (unsigned j = 0, k = 0; j < sPC; j++) {
			for (unsigned i = 0; i < p_nSimP; i++, k++) {				
				io_group3dBasic[k] += mu_est_basic(j+1);
			}
		}
		i_group3dNoisyOrig = i_group3dNoisy;
  }
}


/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeAggregation(
	std::vector<double> &io_im
,	std::vector<double> &io_weight
,	std::vector<bool> &io_mask
,	std::vector<double> const& i_group3d
,	std::vector<unsigned> const& i_index
,	const ImageSize &p_imSize
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls	= p_imSize.nChannels;
	const unsigned width	= p_imSize.width;
	const unsigned wh		= width * p_imSize.height;
	const unsigned sP		= p_params.sizePatch;

	//! Aggregate estimates	
	for (unsigned n = 0; n < p_nSimP; n++) {
		const unsigned ind = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++) {
			const unsigned ij = ind + c * wh;
			for (unsigned p = 0; p < sP; p++) {
				for (unsigned q = 0; q < sP; q++, k++) {
					unsigned ind2 = ij + p * width + q;
					io_im[ind2] += i_group3d[k * p_nSimP + n];
					io_weight[ind2]++;
				}
			}
		}

		//! Apply Paste Trick
		io_mask[ind] = false;

	}
}

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_imReference: image of reference, when the weight if null;
 * io_imResult: will contain the final image;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return : none.
 **/
void computeWeightedAggregation(
	std::vector<double> const& i_imNoisy
,	std::vector<double> &io_imBasic
,	std::vector<double> &io_imFinal
,	std::vector<double> const& i_weight
,	const nlbParams &p_params
,	const ImageSize &p_imSize
){
	for (unsigned c = 0, k = 0; c < p_imSize.nChannels; c++) {

		for (unsigned ij = 0; ij < p_imSize.wh; ij++, k++) {

			//! To avoid weighting problem (particularly near boundaries of the image)
			if (i_weight[k] > 0.) {
				
				io_imFinal[k] /= i_weight[k];
				
			}
			else {
				io_imFinal[k] = io_imBasic[k];
			}
		}
	}
}

