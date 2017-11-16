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
#include "hbe/hbe.h"
#include "Utilities/LibImages.h"

#include "PLE_lib.h"
#include "util.h"
#include "io_pgm.h"

#include <algorithm>

using namespace std;


static void show_help();

/// help on usage of hbe code
static void show_help() {
    cout << "usage: HBE input mask sve_factors PRNU output -muR -sigmaR -tau -gain [options]\n"
         << "1. input: corrupted sve image filename (.pgm)\n"
         << "2. mask image filename (.pgm)\n"
         << "3. sve factors image filename (.pgm)\n"
         << "4. output image filename (.pgm)\n"
         << "Mandatory camera parameters\n"
         << "-muR: camera noise mean value (offset)\n"
         << "-sigmaR: camera noise variance\n"
         << "-tau: shutter speed\n"
         << "-muR: camera gain\n"
         << "Optionnal parameters (default values in parentheses)\n"
         << "-prnu PRNU image filename (.pgm) (default image with 1)\n"
         << "-pSize: patch size (4)\n"
         << "-offset: offset between patches (pSize -1 = 3)\n"
         << "-alpha_H: (1)\n"
         << "-alpha_L: (0.5)\n"
         << "-epsilon_pd: (0.001)\n"
         << std::endl;
}

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 *
 * @author MARC LEBRUN  <marc.lebrun.ik@gmail.com>
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


int main(int argc, char **argv)
{
    //! Check if there is the right call for the algorithm
    if (argc < 8) {

        if(cmdOptionExists(argv, argv+argc, "-help"))
        {
            show_help();
            return 0;
        }
        cout << "usage: HBE input mask sve_factors PRNU output -muR -sigmaR -tau -gain [options]\n"
             << "1. input: corrupted sve image filename (.pgm)\n"
             << "2. mask image filename (.pgm)\n"
             << "3. sve factors image filename (.exr)\n"
             << "4. output image filename (.exr)\n"
             << "Mandatory camera parameters\n"
             << "-muR: camera noise mean value (offset)\n"
             << "-sigmaR: camera noise variance\n"
             << "-tau: shutter speed\n"
             << "-muR: camera gain\n"
             << "Optionnal parameters (default values in parentheses)\n"
             << "-prnu PRNU image filename (.exr) (default image with 1)\n"
             << "-pSize: patch size (4)\n"
             << "-offset: offset between patches (pSize -1 = 3)\n"
             << "-alpha_H: (1)\n"
             << "-alpha_L: (0.5)\n"
             << "-epsilon_pd: (0.001)\n"
             << std::endl;
        return EXIT_FAILURE;
    }

    //! Camera parameters
    double muR = 2046.0; //! camera noise mean value (offset)
    double sigmaR = 30.0;     //! camera noise variance
    double tau = 0.051;  //! shutter speed
    double gain = 0.87;     //! camera gain

    //! Variables initialization
    double alphaH = 1;
    double alphaL = 0.5;
    unsigned pSize = 8;
    double minPixKnown = 0.5;
    unsigned offset = pSize - 1;
    double Nfactor = 1.5;
    double epsilon_pd = 0.1;
    double NfactorPrior = 2.5;
    bool default_prnu = true;
    char * prnu_name = "prnu_tmp.exr";

    if(cmdOptionExists(argv, argv+argc, "-prnu"))
    {
        prnu_name = getCmdOption(argv, argv + argc, "-prnu");
        default_prnu = false;
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

    if(cmdOptionExists(argv, argv+argc, "-muR"))
    {
        muR = (double)atof(getCmdOption(argv, argv + argc, "-muR"));
    } else {
        show_help();
        return EXIT_FAILURE;
    }
    if(cmdOptionExists(argv, argv+argc, "-gain"))
    {
        gain = (double)atof(getCmdOption(argv, argv + argc, "-gain"));
    } else {
        show_help();
        return EXIT_FAILURE;
    }
    if(cmdOptionExists(argv, argv+argc, "-sigmaR"))
    {
        sigmaR = (double)atof(getCmdOption(argv, argv + argc, "-sigmaR"));
    } else {
        show_help();
        return EXIT_FAILURE;
    }
    if(cmdOptionExists(argv, argv+argc, "-tau"))
    {
        tau = (double)atof(getCmdOption(argv, argv + argc, "-tau"));
    } else {
        show_help();
        return EXIT_FAILURE;
    }
    const bool verbose  = true;

    if (verbose) {
        cout << "Parameters: " << endl;
        cout << "Patch size: " << pSize << endl;
        cout << "Patch offset: " << offset << endl;
        cout << "alpha_H: " << alphaH << endl;
        cout << "alpha_L: " << alphaL << endl;
        cout << "Epsilon pd: " << epsilon_pd << endl;
    }

    //! Declarations
    vector<double> imNoisy, imBasic, imFinal, imUmask, imSveFactors;
    ImageSize imSize, imSizeBasic, imSizeUmask, imSizeSve;

    //! Load input image
    if(loadGrayImage_16bits(argv[1], imNoisy, imSize, verbose) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }
    //! Load degradation mask image
    if (strcmp((const char*)(getFileExt(argv[2]).c_str()),"png")==0){
        if(loadImage(argv[2], imUmask, imSizeUmask, verbose) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    } else {
        if(loadGrayImage_16bits(argv[2], imUmask, imSizeUmask, verbose) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    //! Load sve factors image
    if(loadImageExr(argv[3], imSveFactors, imSizeSve, verbose) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }


    if(default_prnu){
        vector<double> data (imSize.whc, 1.);
        saveImageExr(prnu_name, data, imSize);
    }

    //! Load initialization image
    const char * input = argv[1];
    const char * u_mtx = argv[2];
    const char * f_mtx = argv[3];
    const char * prnu_mtx = prnu_name;

    int overlap = offset;
    char * output = "init_tmp.exr";

    int num_orient = 18;
    double epsilon = 0.1;
    double numChannels = 4.0;

    if( PLE_main_routine( 	/* input */ input,
                            /* output */ output,
                            /* u matrix */ u_mtx,
                            /* f matrix */ f_mtx,
                            /* prnu matrix */ prnu_mtx,
                            /* overlap */ overlap,
                            /* patch_size */ pSize,
                            /* num_orientations */ num_orient,//18,
                            /* muR */ muR,
                            /* gain */ gain,
                            /* sigma2R */ sigmaR,
                            /* epsilon for covariance matrix */ epsilon,
                            numChannels) != EXIT_SUCCESS )
        return EXIT_FAILURE;

    if(loadImageExr(output, imBasic, imSizeBasic, verbose) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }



    remove( output );
    if(default_prnu)
        remove( prnu_name );

    //! Normalize image
    for (unsigned k = 0; k < imSize.whc; k++) {
        double factor = gain*tau*imSveFactors[k];
        imNoisy[k] =  imUmask[k]*( imNoisy[k] - muR ) / factor ;
    }

    //! Apply HBE for restoration
    if (verbose) {
        cout << endl << "Restoring the input image usign HBE:" << endl;
    }
    struct timeval start, end;
    gettimeofday(&start, NULL);

    if (runHBE(imNoisy, imBasic, imFinal, imUmask, imSveFactors, imSize,
               verbose, alphaH, alphaL, minPixKnown,
               pSize, offset,Nfactor, NfactorPrior,
               epsilon_pd, muR, sigmaR, gain, tau)!= EXIT_SUCCESS) {
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

    //! Save output image in EXR format
    saveImageExr(argv[4], imFinal, imSize);


    if (verbose) {
        cout << "done." << endl;
    }

    return EXIT_SUCCESS;
}
