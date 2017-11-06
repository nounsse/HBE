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
 * @file utilities.cpp
 * @brief Utilities functions.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "Utilities.h"
#include "LibMatrix.h"

#include <math.h>
#include <omp.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>

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

#ifndef RAND_H
#define RAND_H
extern "C"{
    #include "randmt.h"	//external random number generator
}
#endif

using namespace std;

/**
 * @brief Convenient function to use the sort function provided by the vector library.
 **/
bool comparaisonFirst(
    const pair<double, unsigned> &i_pair1
,	const pair<double, unsigned> &i_pair2
){
    return i_pair1.first < i_pair2.first;
}


NColV centerData(
    std::vector<double> &io_group2d
,	std::vector<double> &i_Umask3d
,	NColV mu_ple
,	const unsigned p_rows
,	const unsigned p_cols
,	double kappa
,	unsigned nSimPreal
,	std::vector<double> &Wout
){
    unsigned NsimP = p_rows;
    unsigned SP2 = p_cols;
    NColV y_rec_cum(SP2);
    y_rec_cum = 0.0;
    NMatrix U_rec_cum(SP2,SP2);
    U_rec_cum = 0.0;

    /* From the index value find the PLE prior */
    NMatrix W(SP2,SP2);

    for (unsigned i = 0; i < nSimPreal; i++) {

        NDiag U(SP2);
        NColV y(SP2);

        unsigned k = 0;
        for (unsigned j = 1; j <= SP2; j++) {
            unsigned ind = (j-1) * NsimP + i;
            /* Recontruct the y vector, which is already centered */
            y(j) = io_group2d[ind];
            /* Reconstruct the U mask for the current patch */
            U(j) = i_Umask3d[ind];

                for (unsigned q = 1; q <= SP2; q++,k++) {
                    W(j,q) = Wout[SP2*SP2*i + k];
                }
        }

        NColV y_rec = W*y;
        y_rec_cum = y_rec_cum + y_rec;
        NMatrix U_rec = W*U;
        U_rec_cum = U_rec_cum + U_rec;
    }

    NColV mu_0(SP2);
    mu_0 = mu_ple;

    NId Id(SP2);

    double kappa_curr = kappa*nSimPreal;
    NColV mu = (kappa_curr*Id + U_rec_cum).i()*(y_rec_cum + kappa_curr*mu_0);

    /* Center data */
    for (unsigned i = 0; i < nSimPreal; i++) {
        for (unsigned j = 0; j < SP2; j++) {
            io_group2d[j * NsimP + i] -= i_Umask3d[j * NsimP + i]*mu(j+1);
        }
    }

    return mu;
}

NColV centerDataForPrior(
    std::vector<double> &io_group3d
,	const unsigned p_rows
,	const unsigned p_cols
){
    NColV muV(p_cols);
    const double inv = 1.0 / (double) p_rows;
    for (unsigned j = 0; j < p_cols; j++) {
        double sum = 0.0;
        for (unsigned i = 0; i < p_rows; i++) {
            sum += io_group3d[j * p_rows + i];
        }

        double mu = sum * inv;

        //! Save mu for output
        muV(j+1) = mu;

        //! Center data
        for (unsigned i = 0; i < p_rows; i++) {
            io_group3d[j * p_rows + i] -= mu;
        }
    }
    return muV;
}


void computeW(
    std::vector<double> &io_group2d
,	std::vector<double> &i_Umask3d
,	std::vector<double> &Wout
,	NSym cov_ple
,	const unsigned p_rows
,	const unsigned p_cols
,	unsigned nSimPreal
,	NDiag Ew
){
    unsigned NsimP = p_rows;
    unsigned SP2 = p_cols;

    NSym Linv(SP2);
    Linv = cov_ple;

    for (unsigned i = 0; i < nSimPreal; i++) {

        NDiag U(SP2);
        NColV y(SP2);

        for (unsigned j = 1; j <= SP2; j++) {
            unsigned ind = (j-1) * NsimP + i;
            /* Recontruct the y vector, which is already centered */
            y(j) = io_group2d[ind];
            /* Reconstruct the U mask for the current patch */
            U(j) = i_Umask3d[ind];
        }

        NMatrix LinvU = Linv*U;
        NMatrix aux = U*LinvU + Ew;
        NSym ULinvUtEw(SP2);
        ULinvUtEw << aux;

        NSym Letoile = ULinvUtEw.i();
        NMatrix W = LinvU*Letoile;

        unsigned k = 0;
        for (unsigned p = 1; p <= SP2; p++)
            for (unsigned q = 1; q <= SP2; q++,k++)
                Wout[ SP2*SP2*i + k ] = W(p,q);
    }

}

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
){
    return (i_value < i_min ? i_min : (i_value > i_max ? i_max : i_value));
}



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
){
    if (i_n == 1) {
        o_a = 1;
        o_b = 1;
        return;
    }

    o_b = 2;
    while (i_n % o_b > 0) {
        o_b++;
    }
    o_a = i_n / o_b;

    if (o_b > o_a) {
        o_a = o_b;
        o_b = i_n / o_a;
    }
}

NSym *read_covs( int num_models, int patch_size ) {

    int data_size = patch_size*patch_size;
    NSym temp( data_size );
    NSym * covs = new NSym [ num_models ];
    for( int model = 0; model < num_models; model++ )
        covs[ model ] = temp;

    FILE *fp = fopen("covs.txt","r");
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

NColV *read_mus( int num_models, int patch_size ) {

    int data_size = patch_size * patch_size;
    NColV * mus = new NColV [num_models];
    NColV mu( data_size );
    mu = 0;
    for( int i = 0; i < num_models; i++ )
        mus[i] = mu;

    FILE *fp = fopen("mus.txt","r");
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

void write_covs( NSym * covs, int num_models, int patch_size ) {

    int data_sz = patch_size*patch_size;
    FILE *fp = fopen("covs2.txt","w");

    for (int c=0; c < num_models; c++) {

        for (int i=1; i<= data_sz; i++)
            for (int j=1; j<= data_sz; j++)
                fprintf(fp,"%f ", covs[c](i,j));
    }

    fclose(fp);
}

void write_mus( NColV * mus, int num_models, int patch_size ) {

    int data_sz = patch_size*patch_size;
    FILE *fp = fopen("mus2.txt","w");

    for (int c=0; c < num_models; c++) {

        for (int i=1; i<= data_sz; i++)
                fprintf(fp,"%f ", mus[c](i));
    }

    fclose(fp);
}

void write_mu( const char *fname, NColV mu, int patch_size ) {

    int data_sz = patch_size*patch_size;
    FILE *fp = fopen(fname,"w");

    for (int i=1; i<= data_sz; i++)
                fprintf(fp,"%f ", mu(i));

    fclose(fp);
}

void write_cov( const char *fname, NSym covs, int patch_size ) {

    int data_sz = patch_size*patch_size;
    FILE *fp = fopen(fname,"w");

    for (int i=1; i<= data_sz; i++)
        for (int j=1; j<= data_sz; j++)
            fprintf(fp,"%lf ", covs(i,j));

    fclose(fp);
}

void write_matrix( const char *fname, NMatrix mtx, int patch_size ) {

    int data_sz = patch_size*patch_size;
    FILE *fp = fopen(fname,"w");
    //FILE *fp = fopen("mtx.txt","w");

    for (int i=1; i<= data_sz; i++)
        for (int j=1; j<= data_sz; j++)
            fprintf(fp,"%lf ", mtx(i,j));

    fclose(fp);
}


void write_mask( const char *fname, vector<bool> mask, int size ) {

    FILE *fp = fopen(fname,"w");

    for (int i=0; i< size; i++)
            if (mask[i])
                fprintf(fp,"%lf ", 1.0);
            else
                fprintf(fp,"%lf ", 0.0);

    fclose(fp);
}

void write_patch( const char *fname, std::vector<double> const& image, unsigned im_size, unsigned sP, unsigned ij ) {

    FILE *fp = fopen(fname,"w");

    for (unsigned p = 0; p < sP; p++)
        for (unsigned q = 0; q < sP; q++)
            fprintf(fp,"%lf ", image[ij + p * im_size + q]);

    fclose(fp);

}

void write_group3d( const char *fname, std::vector<std::vector<double> > &i_group3d, unsigned sP, unsigned nSimP ) {

    FILE *fp = fopen(fname,"w");

    unsigned total = sP*sP*nSimP;

    for (unsigned p = 0; p < total; p++)
        fprintf(fp,"%lf ", i_group3d[0][p]);

    fclose(fp);
}

void write_group2d( const char *fname, std::vector<double> &i_group2d, unsigned sP, unsigned nSimP ) {

    FILE *fp = fopen(fname,"w");

    unsigned total = sP*sP*nSimP;

    for (unsigned p = 0; p < total; p++)
        fprintf(fp,"%lf ", i_group2d[p]);

    fclose(fp);
}

unsigned pixels_to_fillMask( vector<bool> mask, unsigned size) {

    unsigned toFill = 0;

    for (unsigned ij = 0; ij < size; ij++)
        if (mask[ij])
            toFill++;

    return toFill;

}

unsigned pixels_to_fill( vector<bool> mask, std::vector<double> &i_Umask, unsigned size) {

    unsigned toFill = 0;

    for (unsigned ij = 0; ij < size; ij++)
        if ((i_Umask[ij]==0) && mask[ij] )
            toFill++;

    return toFill;

}

















