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
 * @file LibMatrix.cpp
 * @brief Tools for matrix manipulation, based on ccmath functions
 *        by Daniel A. Atkinson.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibMatrix.h"
#include "../Utilities/Utilities.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
 
#ifndef SHORT_NEWMAT
#define SHORT_NEWMAT
#include "newmatap.h"  // NEWMAT 
typedef NEWMAT::Matrix NMatrix;
typedef NEWMAT::SymmetricMatrix NSym;
typedef NEWMAT::DiagonalMatrix NDiag;
typedef NEWMAT::IdentityMatrix NId;
typedef NEWMAT::ColumnVector NColV;
typedef NEWMAT::RowVector NRowV;
#endif
typedef NEWMAT::LowerTriangularMatrix NLowT;
typedef NEWMAT::UpperTriangularMatrix NUpT;

 using namespace std;


NSym covarianceMatrixForPrior(
	std::vector<double> &io_group2d
,	const unsigned NsimP
,	const unsigned SP2
,	double epsilon_pd
){

	NMatrix cov(SP2,SP2);
	cov = 0.0;

	for (unsigned i = 0; i < NsimP; i++) {
		
		NColV y(SP2);
		for (unsigned j = 0; j < SP2; j++) {
			//! data is already centered
			y(j+1) = io_group2d[j * NsimP + i]; 				
		}		
		
		cov = cov + y*y.t();
	}
	cov = cov / (NsimP-1);
	
	NSym covariance;
	NId Id(SP2);
	covariance << cov + Id*epsilon_pd;
	
	return covariance;
}


NSym covarianceMatrix(
	std::vector<double> &io_group2d
,	std::vector<double> &i_Umask3d
,	NColV mu_ple
,	NSym cov_ple
,	const unsigned NsimP
,	const unsigned SP2
,	double kappa
,	double nu
,	NColV mu // mu estimated in centerData function
,	double beta
,	unsigned nSimPreal
,	double epsilon_pd
,	NDiag Ew
){
	NMatrix y_rec_cum(SP2,SP2);
	NMatrix Linv(SP2,SP2);
	y_rec_cum = 0.0;
		
	for (unsigned i = 0; i < nSimPreal; i++) {	
			NColV y(SP2);
			for (unsigned j = 1; j <= SP2; j++) {
				unsigned ind = (j-1) * NsimP + i;
				y(j) = io_group2d[ind];	//! values already centered				
			}										
			y_rec_cum = y_rec_cum + y*y.t();	
	}	

	double nu_curr = SP2 + nu*nSimPreal;	
	double kappa_curr = kappa*nSimPreal;
	
	NColV mu_0(SP2);
	mu_0 = mu_ple;
	NSym Sigma_0 = cov_ple;
		
	NMatrix L_hatInv = (Sigma_0*nu_curr + kappa_curr*(mu - mu_0)*
					   (mu - mu_0).t() + y_rec_cum) / 
					   (nu_curr + nSimPreal - SP2);			
	NId Id(SP2);
	Linv = (L_hatInv + L_hatInv.t())/2.0 + Id*epsilon_pd;												
	
	NSym cov_out(SP2);
	cov_out << Linv;
	
	return cov_out;
}



