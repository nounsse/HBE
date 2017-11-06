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
#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <vector>
#include "LibImages.h"

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


/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
    const std::pair<double, unsigned> &i_pair1
,	const std::pair<double, unsigned> &i_pair2
);

/**
 * @brief Clip a value between min and max
 *
 * @param i_value: value to clip;
 * @param i_min: minimum value;
 * @param i_max: maximum value.
 *
 * @return value clipped between [min, max].
 **/
double clip(
    const double i_value
,	const double i_min
,	const double i_max
);

/**
 * @brief Obtain and substract the baricenter of io_group3d.
 *
 * @param io_group3d(p_rows x p_cols) : data to center;
 * @param o_baricenter(p_cols): will contain the baricenter of io_group3d;
 * @param p_rows, p_cols: size of io_group3d.
 *
 * @return none.
 **/
void computeW(
    std::vector<double> &io_group2d
,	std::vector<double> &i_Umask3d
,	std::vector<double> &Wout
,	NSym cov_ple
,	const unsigned p_rows
,	const unsigned p_cols
,	unsigned nSimPreal
,	NDiag Ew
);

NColV centerData(
    std::vector<double> &io_group2d
,	std::vector<double> &i_Umask3d
,	NColV mu_ple
,	const unsigned p_rows
,	const unsigned p_cols
,	double kappa
,	unsigned nSimPreal
,	std::vector<double> &Wout
);

NColV centerDataForPrior(
    std::vector<double> &io_group3d
,	const unsigned p_rows
,	const unsigned p_cols
);



/**
 * @brief Determine a and b such that : n = a * b, with a and b as greatest as possible
 *
 * @param i_n : number to decompose;
 * @param o_a : will contain a;
 * @param o_b : will contain b.
 *
 * @return none.
 **/
void determineFactor(
    const unsigned i_n
,   unsigned &o_a
,   unsigned &o_b
);


NSym * read_covs( int num_models, int patch_size );

NColV * read_mus( int num_models, int patch_size );

void write_covs( NSym * covs, int num_models, int patch_size );

void write_mus( NColV * mus, int num_models, int patch_size );
void write_mu( const char *fname, NColV mu, int patch_size );
void write_cov( const char *fname, NSym covs, int patch_size );
void write_matrix( const char *fname, NMatrix mtx, int patch_size );

void write_mask( const char *fname, vector<bool> mask, int size );
unsigned pixels_to_fillMask( vector<bool> mask, unsigned size);
unsigned pixels_to_fill( vector<bool> mask, std::vector<double> &i_Umask, unsigned size);
void write_group3d( const char *fname, std::vector<std::vector<double> > &i_group3d, unsigned sP, unsigned nSimP );
void write_group2d( const char *fname, std::vector<double> &i_group2d, unsigned sP, unsigned nSimP );
void write_patch( const char *fname, std::vector<double> const& image, unsigned im_size, unsigned sP, unsigned ij );

#endif // UTILITIES_H_INCLUDED
