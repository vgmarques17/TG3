#include "mylib.h"
#include <stdio.h>
#include <stdlib.h>

void showchar(char ch){
	
	printf("Your character was %c", ch);
	return;
}


float *getMatrixCol(int n_lines,int n_cols, float matrix[n_lines][n_cols], int col){
	
	if(col>=n_cols){
		printf("Erro! Matriz possui apenas %d colunas. \n",n_cols);
		return NULL;
	}
	
	float* vector = malloc(n_lines);
	
	for(int i=0;i<n_lines;i++){
		vector[i] = matrix[i][col];
	}
	return vector;
}
