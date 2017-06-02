#ifndef _MYLIB_H_
#define _MYLIB_H_

#include "mylib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#define PI 3.14159265 
#define QUASI_ZERO 1e-15

typedef struct {
	float *amplitude;// Amplitude do seno modulado para cada canal
	float *phase_rad;// Fase do...
	float *offset;// offset do ...
	float *alpha;// 
	float *beta;//
	float *omegat;//
	float *y;//
}fitSine;

/*Intern functions*/



/*Extern functions*/
extern float * getMatrixCol(int n_lines,int n_cols, float matrix[n_lines][n_cols], int col); //returns a column of file as an array

extern void showTable_daq(float* table, int rows, int col); // Plots a table on the screen "1D" pointers

extern void showTable(float** table,int Nrows, int Ncol); //Plots a table on the screen, for "2D" pointers

extern void saveData(float* Data, int rows, int col); // Save a table to a file, ';' separated, no header

extern int pointer2gsl(gsl_matrix * gsl_M, float* M, int m, int n); //converter dados para matriz GSL

extern int gsl2pointer(float* M, gsl_matrix* gsl_M);//converter matriz GSL para organizacao com pointers

int gsl2pointer_vector(float* M, gsl_vector* gsl_M); //converte vetor gsl para pointers

extern int transposeMatrix(float **usr_matrix, int Nrow,int Ncol, float **T_matrix); // transpose using pointers

extern int productMatrix(float** A, int rA, int cA,float ** B,int rB, int cB,float** AB);// matrix product using pointers

extern int gsl_productMatrix(gsl_matrix *A, gsl_matrix *B, gsl_matrix * AB);// matrix product using GSL

extern int pinv(gsl_matrix *usr_matrix,gsl_matrix *P_inv); //pseudo inverse algorithm using SVD

extern fitSine sineRegression_lms(gsl_matrix* Data, float f0,float samp_freq); // fixed frequency sine regression, talvez n precise de k0




#endif
