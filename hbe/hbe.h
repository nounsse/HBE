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
#ifndef HBE_H_INCLUDED
#define HBE_H_INCLUDED

#include<vector>

#include "../Utilities/Utilities.h"

/**
 * @brief Structures of parameters dedicated to HBE process
 *
 * @param sigma: value of noise;
 * @param sizePatch: size of patches (sizePatch x sizePatch);
 * @param nSimilarPatches: minimum number of similar patches wanted;
 * @param sizeSearchWindow: size of the search window around the reference patch;
 * @param boundary: must be > sizeSearchWindow. Boundary kept around sub-images when the image is
 *      subdivided for parallelization;
 * @param offSet: step between two similar patches;
 * @param useHomogeneousArea: if true, use the homogeneous area trick;
 * @param gamma: threshold to detect homogeneous area;
 * @param beta: parameter used to estimate the covariance matrix;
 * @param tau: parameter used to determine similar patches;
 * @param isFirstStep: true if the first step of the algorithm is wanted;
 * @param doPasteBoost: if true, patches near denoised similar patches will not be used as reference
 *		patches;
 * @param verbose: if true, print informations.
 **/
struct nlbParams
{
	double sigma;
	unsigned sizePatch;
	unsigned nSimilarPatches;
	unsigned sizeSearchWindow;
	unsigned boundary;
	unsigned offSet;
	bool useHomogeneousArea;
	double gamma;
	double beta;
	double tau;
	bool isFirstStep;
	bool doPasteBoost;
	bool verbose;
};

/**
 * @brief Structure containing usefull matrices for the Bayes estimations.
 *
 * @param group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 * @param baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 * @param covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 * @param covMatTmp: allocated memory. Used to process the Bayes estimate;
 * @param tmpMat: allocated memory. Used to process the Bayes estimate.
 **/
struct matParams
{
	std::vector<double> group3dTranspose;
	std::vector<double> baricenter;
	std::vector<double> covMat;
	std::vector<double> covMatTmp;
	std::vector<double> tmpMat;
};

/**
 * @brief Initialize Parameters of the HBE algorithm.
 *
 * @param o_paramStep1 : will contain the nlbParams for the first step of the algorithm;
 * @param o_paramStep2 : will contain the nlbParams for the second step of the algorithm;
 * @param p_sigma : standard deviation of the noise;
 * @param p_imSize: size of the image;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_verbose : if true, print some informations.
 *
 * @return none.
 **/
void initializeParameters(
	nlbParams &o_paramStep2
,   const double p_sigma
,	const ImageSize &p_imSize
,   const bool p_verbose
,	unsigned pSize
,	unsigned offset
);

/**
 * @brief Main function to process the whole HBE algorithm.
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
);

/**
 * @brief Generic step of the HBE denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_imSize: size of i_imNoisy;
 * @param p_params: parameters of the method, contains:
 *			- sigma: standard deviation of the noise;
 *			- sizePatch: size of patches (sizePatch x sizePatch);
 *			- nSimilarPatches: number of similar patches;
 *			- sizeSearchWindow: size of the neighbourhood searching window;
 *			- useHomogeneousArea: if true, the trick of using homogeneous area will be used;
 *			- gamma: parameter used to determine if we are in an homogeneous area;
 *			- maxAvoid: parameter used to stop the paste trick;
 *			- beta: parameter used during the estimate of the denoised patch;
 *			- coefBaricenter: parameter to determine if the covariance matrix inversion is correct;
 *			- isFirstStep: true if it's the first step of the algorithm which is needed;
 *			- verbose: if true, print some informations, do nothing otherwise.
 *
 * @return none.
 **/
void processHBE(
	std::vector<double> const& i_imNoisy
,	std::vector<double> &io_imBasic
,	std::vector<double> &i_Umask
,	std::vector<double> &o_imFinal
,	const ImageSize &p_imSize
,	nlbParams &p_params
,	double alphaH
,	double alphaL
,	double minPixKnown
,	double Nfactor
,	double NfactorPrior
,	double epsilon_pd
,	double varSigma
);


/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param i_imBasic: contains the basic estimation;
 * @param o_group3dNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_group3dBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_imSize: size of images;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned findSimilarPatches(
	std::vector<double> const& i_imNoisy
,	std::vector<double> const& i_imBasic
,	std::vector<double> const& i_Umask
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
,	unsigned *nSimPprior
);


/**
 * @brief Compute the Bayes estimation.
 *
 * @param i_group3dNoisy: contains all similar patches in the noisy image;
 * @param io_group3dBasic: contains all similar patches in the basic image. Will contain estimates
 *			for all similar patches;
 * @param i_mat: contains :
 *		- group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_imSize: size of the image;
 * @param p_params: see processStep2 for more explanations;
 * @param p_nSimP: number of similar patches.
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
);


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
);

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_imReference: image of reference, when the weight if null;
 * io_imResult: will contain the final image;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return none.
 **/
void computeWeightedAggregation(
	std::vector<double> const& i_imNoisy
,	std::vector<double> &io_imBasic
,	std::vector<double> &io_imFinal
,	std::vector<double> const& i_weight
,	const nlbParams &p_params
,	const ImageSize &p_imSize
);


#endif // HBE_H_INCLUDED
