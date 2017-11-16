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
 * @file PLE_denoise.cpp
 * @brief command-line handler for PLE image denoising
 * @author Yi-Qing WANG <yiqing.wang@polytechnique.edu>
 */


#ifndef STD_H
#define STD_H
#include <iostream>
#include <iomanip>
#include <cstdio>
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif

#ifndef PLE_H
#define PLE_H
#include "PLE_lib.h" 		//routines dedicated to this implementation
#endif

using namespace std;

#if 0
int main( int argc, char *const *argv ){


	if( !(argc == 12) ){

		cerr << endl << endl << endl;
		cerr << "\t\t#######################################################################################################" << endl << endl;
		cerr << "\t\tPLE initialization for HDR image generation" << endl << endl;
		cerr << "\t\t" << endl << endl;
		cerr << "\t\t#######################################################################################################" << endl << endl;
		cerr << "\t Usage : please specify the parameters for " << argv[0] << endl << endl; 
		cerr << "\t\t 1. corrupted sve image filename (.pgm)" << endl;
		cerr << "\t\t 2. mask image filename (.pgm)" << endl;
		cerr << "\t\t 3. sve factors image filename (.pgm)" << endl;				
		cerr << "\t\t 4. PRNU image filename (.pgm)" << endl;
		cerr << "\t\t " << endl;
		cerr << "\t Camera parameters:" << endl;
		cerr << "\t\t " << endl;
		cerr << "\t\t 5. noise offset" << endl;
		cerr << "\t\t 6. noise variance" << endl;
		cerr << "\t\t 7. camera gain" << endl;
		cerr << "\t\t 8. shutter speed" << endl;
		cerr << "\t\t 9. patch size" << endl;
		cerr << "\t\t 10. patch overlap" << endl;							
		cerr << "\t\t 11. output image (EXR format)" << endl;																			
		cerr << "\t\t " << endl;
		
        	return EXIT_FAILURE;
	}		


	const char * input = argv[1];
	const char * u_mtx = argv[2];
	const char * f_mtx = argv[3];
	const char * prnu_mtx = argv[4];
	double muR = atof(argv[5]);
	double sigma2R = atof(argv[6]);
	double g = atof(argv[7]);	
	double tau = atof(argv[8]);

	int pSize = atoi(argv[9]);
	int overlap = atoi(argv[10]);				
	const char * output = argv[11];
	
	int num_orient = 18;
	double epsilon = 0.1;	
	double numChannels = 4.0;

	return PLE_main_routine( 	/* input */ input,
					/* output */ output,
					/* u matrix */ u_mtx,
					/* f matrix */ f_mtx,
					/* prnu matrix */ prnu_mtx,
					/* overlap */ overlap, 
					/* patch_size */ pSize,
					/* num_orientations */ num_orient,//18,
					/* muR */ muR,
					/* gain */ g,
					/* sigma2R */ sigma2R,
					/* epsilon for covariance matrix */ epsilon,
					numChannels);
}
#endif
