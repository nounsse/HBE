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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <string>
#include <cstring>
#include <sstream>

#include <sys/time.h>

#include "Utilities/Utilities.h"
#include "PLE_lib.h"
#include "Utilities/LibImages.h"
#include "hbe.h"

#include <algorithm>

using namespace std;


static void show_help();

/// help on usage of hbe code
static void show_help() {
    std::cerr<<"\nHBE.\n"
            << "Usage: "
            << " HBE input mask initImage output [options]\n\n"
            << "Options (default values in parentheses)\n"
            << "-sigma: Degradation parameter (30)\n"
            << "-pSize: patch size (4)\n"
            << "-offset: offset between patches (pSize -1 = 3)\n"
            << "-alpha_H: (1)\n"
            << "-alpha_L: (0.5)\n"
            << "-epsilon_pd: (0.001)\n"
            << "-noise: (0 : no noise, 1 : noise, default : 0)"
            << std::endl;
}

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 *
 **/


/**
 *
 */
/**
* @brief Find the command option named option
*/
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

/**
 *
 */
/**
* @brief Get file exension of file name
*
* @param File name
* @return File extension
*/
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

/**
 *
 */
/**
* @brief Check for input parameter
*
* @param beginning of input command chain
* @param end of input command chain
* @return whether the parameter exists or not
*/
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

// Computes random number between min and max
// between 0 and 1 if min == max
double getRandomNumber( double min, double max ){

    double min_v = std::min( min, max );
    double max_v = std::max( min, max );

    if( min_v == max_v ){
        min_v = 0;
        max_v = 1;
    }

    double r = ((double) rand() / (RAND_MAX));
    return ( r )*( max_v - min_v ) + min_v;
}

int main(int argc, char **argv)
{
    //! Check if there is the right call for the algorithm
    if (argc < 3) {

        if(cmdOptionExists(argv, argv+argc, "-help"))
        {
            show_help();
            return 0;
        }

        cout << "usage: HBE input output [options]\n\n"
             << "Options (default values in parentheses)\n"
             << "-sigma: Degradation parameter (30)\n"
             << "-pSize: patch size (4)\n"
             << "-offset: offset between patches (pSize -1 = 3)\n"
             << "-alpha_H: (1)\n"
             << "-alpha_L: (0.5)\n"
             << "-epsilon_pd: (0.001)\n"
             << "-degradation: generate a degradation mask of a given percentage\n"
             << "-mask: degradation mask taken from a file filename.png\n"
             << "-noise: (0 : no noise, 1 : noise, default : 0)"
             << std::endl;
        return EXIT_FAILURE;
    }

    //! Variables initialization
    double sigma = 30;
    double alphaH = 1;
    double alphaL = 0.5;
    unsigned pSize = 4;
    double minPixKnown = 0.5;
    unsigned offset = pSize - 1;
    double Nfactor = 1.5;
    double epsilon_pd = 0.001;
    double NfactorPrior = 2.5;
    unsigned doNoise = 0;
    char * mask_filename;
    bool use_mask_file =false;
    bool use_degradation_mask = false;
    int degradation = 0;
    char * degradated_img = argv[1];

    if(cmdOptionExists(argv, argv+argc, "-sigma"))
    {
        sigma = (double)atof(getCmdOption(argv, argv + argc, "-sigma"));
    }
    if(cmdOptionExists(argv, argv+argc, "-pSize"))
    {
        pSize = atoi(getCmdOption(argv, argv + argc, "-pSize"));
    }
    if(cmdOptionExists(argv, argv+argc, "-offset"))
    {
        offset = atoi(getCmdOption(argv, argv + argc, "-offset"));
    } else {
        offset = pSize - 1;
    }
    if(cmdOptionExists(argv, argv+argc, "-alpha_H"))
    {
        alphaH = (double)atof(getCmdOption(argv, argv + argc, "-alpha_H"));
    }
    if(cmdOptionExists(argv, argv+argc, "-alpha_L"))
    {
        alphaL = (double)atof(getCmdOption(argv, argv + argc, "-alpha_L"));
    }
    if(cmdOptionExists(argv, argv+argc, "-epsilon_pd"))
    {
        epsilon_pd = (double)atof(getCmdOption(argv, argv + argc, "-epsilon_pd"));
    }
    if(cmdOptionExists(argv, argv+argc, "-noise"))
    {
        doNoise = atoi(getCmdOption(argv, argv + argc, "-noise"));
    }

    if(cmdOptionExists(argv, argv+argc, "-mask"))
    {
        use_mask_file = true;
        mask_filename = getCmdOption(argv, argv + argc, "-mask");
    }

    if(cmdOptionExists(argv, argv+argc, "-degradation"))
    {
        degradation = atoi(getCmdOption(argv, argv + argc, "-degradation"));
    }

    double varSigma = sigma*sigma;

    const bool verbose  = true;

    if (verbose) {
        cout << "Parameters: " << endl;
        cout << "Sigma: " << sigma << endl;
        cout << "Patch size: " << pSize << endl;
        cout << "Patch offset: " << offset << endl;
        cout << "alpha_H: " << alphaH << endl;
        cout << "alpha_L: " << alphaL << endl;
        cout << "Epsilon pd: " << epsilon_pd << endl;
        cout << "Noise: " << doNoise << endl;
    }

    //! Declarations
    vector<double> im, imNoisy, imBasic, imFinal, imUmask;
    ImageSize imSize, imSizeBasic, imSizeUmask;

    //! Load input image
    if (strcmp((const char*)(getFileExt(argv[1]).c_str()),"pgm")==0)
    {
        if(loadGrayImage_16bits(argv[1], im, imSize, verbose) != EXIT_SUCCESS)
            return EXIT_FAILURE;
    } else if (strcmp((const char*)(getFileExt(argv[1]).c_str()),"png")==0){
        if(loadImage(argv[1], im, imSize, verbose) != EXIT_SUCCESS)
            return EXIT_FAILURE;
    } else {
        return EXIT_FAILURE;
    }

    if( use_mask_file ){

        //! Load degradation mask image
        if (strcmp((const char*)(getFileExt(mask_filename).c_str()),"png")==0)
        {
            if(loadImage(mask_filename, imUmask, imSizeUmask, verbose) != EXIT_SUCCESS)
                return EXIT_FAILURE;
        }
        else {
            if(loadGrayImage_16bits(mask_filename, imUmask, imSizeUmask, verbose) != EXIT_SUCCESS)
                return EXIT_FAILURE;
        }
        use_degradation_mask = true;
    } else {

        //! Generate degradation mask
        imSizeUmask.width = imSize.width;
        imSizeUmask.height = imSize.height;
        imSizeUmask.nChannels = 1;
        imSizeUmask.wh =  imSizeUmask.width*imSizeUmask.height;
        imSizeUmask.whc =  imSizeUmask.width*imSizeUmask.height*imSizeUmask.nChannels;

        imUmask.clear();
        imUmask.resize(imSize.wh, 0);

        if (degradation > 0) {
            int count = 0;
            int max_number = imSizeUmask.wh*(float)degradation/100.;

            srand(time(NULL)); // initialisation de rand

            while ( count < max_number ){
                int rand ;
                do {
                    rand = getRandomNumber( 0, imSizeUmask.wh );
                } while( imUmask[rand] != 0 );
                imUmask[rand] = 128;
                count++;
            }

            mask_filename = "tmp_mask_file.png";
            saveImage(mask_filename, imUmask, imSizeUmask, 0, 255);
            use_degradation_mask = true;
        }
    }

    //! Add noise
    if ( doNoise == 1 ) {
        addNoise(im, imNoisy, sigma, verbose);
    }
    else
        imNoisy = im;

    vector<double> degradated = imNoisy;
    for(unsigned int i = 0 ; i < imUmask.size() ; i ++ )
        if( imUmask[i] > 0 )
            degradated[i] = imUmask[i];

    if( doNoise == 1 || use_degradation_mask ){
        degradated_img = "degradated.png";
        saveImage(degradated_img, degradated, imSize, 0, 255);
    }

    //! Load initialization image

    if ( use_degradation_mask ) {
        const char * input = degradated_img;
        const char * u_mtx = mask_filename;
        int overlap = offset;
        char * output = "tmp_init_img.png";

        int num_orient = 18;
        double epsilon = 30.0;
        double numChannels = 1.0;

        PLE_main_routine( 	/* input */ input,
                            /* output */ output,
                            /* u matrix */ u_mtx,
                            /* polluting_sigma */ sigma,
                            /* overlap */ overlap,
                            /* patch_size */ pSize,
                            /* num_orientations */ num_orient,
                            /* epsilon for covariance matrix */ epsilon,
                            numChannels);

        if(loadImage(output, imBasic, imSizeBasic, verbose) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }

        remove( output );
        if(degradation > 0)
            remove( mask_filename );

    } else {
        imBasic = imNoisy;
        imSizeBasic = imSize;
    }

    //! Apply HBE for restoration
    if (verbose) {
        cout << endl << "Restoring the input image using HBE:" << endl;
    }
    struct timeval start, end;
    gettimeofday(&start, NULL);

    if (runHBE(imNoisy, imBasic, imFinal, imUmask, imSize,
               sigma, verbose, alphaH, alphaL, minPixKnown,
               pSize, offset,Nfactor, NfactorPrior,
               epsilon_pd, varSigma)!= EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    if (verbose) {
        gettimeofday(&end, NULL);
        double elapsedTime = (end.tv_sec  - start.tv_sec) +
                (end.tv_usec - start.tv_usec) / 1.e6;
        std::cout<< "time elapsed : "<< elapsedTime << " seconds "<< std::endl;
        std::cout<< "***************************"<< std::endl<< std::endl<< std::endl;

        cout << endl;
    }

    std::string output_filename (argv[2]);

    //! Save output image
    if (strcmp((const char*)(getFileExt(argv[2]).c_str()),"png")!=0)
        output_filename.append(".png");

    saveImage(output_filename.c_str(), imFinal, imSize, 0, 255);

    if (verbose) {
        cout << "done." << endl;
    }

    return EXIT_SUCCESS;
}
