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
#ifndef LIB_MATRIX_H_INCLUDED
#define LIB_MATRIX_H_INCLUDED

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

#include<vector>

/**
 * @brief Compute the covariance matrix.
 *
 * @param i_patches: set of patches of size (nb x N);
 * @param o_covMat: will contain patches' * patches;
 * @param p_N : size of a patch;
 * @param p_nb: number of similar patches in the set of patches.
 *
 * @return none.
 **/
NSym covarianceMatrixForPrior(
	std::vector<double> &io_group2d
,	const unsigned NsimP
,	const unsigned SP2
,	double epsilon_pd
);

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
);



#endif // LIB_MATRIX_H_INCLUDED
