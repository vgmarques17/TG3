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
#include <libconfig.h>

#define PI 3.14159265 
#define QUASI_ZERO 1e-15
#define TRUE 1
#define FALSE 0

typedef struct {
	gsl_vector *amplitude;// Amplitude do seno modulado para cada canal
	gsl_vector *phase_rad;// Fase do...
	gsl_vector *offset;// offset do ...
	gsl_vector *alpha;// 
	gsl_vector *beta;//
	gsl_vector *omegat;//
	gsl_matrix *y;//
}fitSine;

/*Intern functions*/



/*Extern functions*/


extern void showTable_1Df(float* table, int rows, int col); // Plots a table on the screen "1D" pointers

extern void showTable_2Df(float** table,int Nrows, int Ncol); //Plots a table on the screen, for "2D" pointers

extern void saveData_1Df(float* Data, int rows, int col); // Save a table to a file, ';' separated, no header

extern int pointer2gsl_matrix(gsl_matrix * gsl_M, float* M, int m, int n); //converter dados para matriz GSL

extern int gsl_matrix2pointer(float* M, gsl_matrix* gsl_M);//converter matriz GSL para organizacao com pointers

extern int pointer2gsl_array(gsl_vector * gsl_vector, float* M); //converter dados para vetor GSL

extern int gsl_array2pointer(float* M, gsl_vector* gsl_M); //converte vetor gsl para pointers

extern int transposeMatrix(float **usr_matrix, int Nrow,int Ncol, float **T_matrix); // transpose using pointers

extern int productMatrix(float** A, int rA, int cA,float ** B,int rB, int cB,float** AB);// matrix product using pointers

extern int gsl_productMatrix(gsl_matrix *A, gsl_matrix *B, gsl_matrix * AB);// matrix product using GSL

extern int pinv(gsl_matrix *usr_matrix,gsl_matrix *P_inv); //pseudo inverse algorithm using SVD

extern fitSine sineRegression_lms(gsl_matrix* Data, float f0,float samp_freq); // fixed frequency sine regression

extern gsl_matrix * reorganizaDados(gsl_matrix* Data,int samples_cicle, int channel, int phase); //pega somente as partes com sinal nas matrizes

extern gsl_matrix * readFile_gsl(const char *filename, int m, int n, int header); //returns a gsl_matrix read from file (maybe can be done directly), 
																				  //no momento soh ignora o header, ideal seria ler tudo

extern gsl_matrix * makeImpedanceTable(gsl_matrix* Data, int samples_cicle, gsl_matrix * impedance_matrix, float f0, float fs); //makes impedance file

extern int saveFile_gsl(config_t configuration, gsl_matrix * Data, int fileType);

#endif
