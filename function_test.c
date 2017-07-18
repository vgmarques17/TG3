#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>
#include <libconfig.h>
#include "../pmd.h"
#include "../usb-1608FS-Plus.h"
#include "Library/mylib.h"

#define MAX_COUNT     (0xffff)
#define FALSE 0
#define TRUE 1

int main(int argc, char **argv){
	
	int L = 1000, factor = 10;
	float f0 = 1000,fs = 5000, phi1 = PI/2, phi2 = PI/4;
	gsl_vector  *seno1, *seno2, *seno3, *seno4;
	gsl_matrix * data, *impedance;

	config_t cfg_params;
	config_setting_t *setting;
	char cfg_file[100];
	
	
	seno1 = gsl_vector_alloc(L);
	seno2 = gsl_vector_alloc(L);
	seno3 = gsl_vector_alloc(L);
	seno4 = gsl_vector_alloc(L);
	
	config_init (&cfg_params);	
	printf("Enter .cfg file name/path:\n");
	scanf("%123s",cfg_file);
	strcat(cfg_file,".cfg");
	
	/* Read the file. If there is an error, report it and exit. */
	if(! config_read_file(&cfg_params, cfg_file))
	{
	fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg_params),
			config_error_line(&cfg_params), config_error_text(&cfg_params));
	config_destroy(&cfg_params);
	return(EXIT_FAILURE);
	}
	
	
	
	for(int i=0;i<L;i++){
		gsl_vector_set(seno1, i,1*sin(2*PI*f0/fs*i+phi2));
		gsl_vector_set(seno2,i, 2*sin(2*PI*f0/fs*i+phi1)+1);	
		gsl_vector_set(seno3, i,3*sin(2*PI*f0/fs*i+phi1+phi2));
		gsl_vector_set(seno4,i, 4*sin(2*PI*f0/fs*i+phi1+phi2)+1);		
		
		}
		
	data = gsl_matrix_alloc(L,4);
	gsl_matrix_set_col(data,0,seno1);
	gsl_matrix_set_col(data,1,seno2);
	gsl_matrix_set_col(data,2,seno3);
	gsl_matrix_set_col(data,3,seno4);
	
	impedance = gsl_matrix_alloc(factor,12);
	
	makeImpedanceTable(data, L/factor, impedance, f0, fs);
	
	for(int i=0;i<factor;i++){
		for(int j = 0;j<12;j++){
			printf("%f | ",gsl_matrix_get(impedance,i,j));
		}
		printf("\n");
	}

	saveFile_gsl(cfg_params,impedance, 1);
	
	return 0;

}

// gcc -g -Wall -I. -o function_test function_test.c -L. Library/mylib.o -lmccusb -L/usr/local/lib -lhidapi-libusb -lusb-1.0 -lm -lgsl -lgslcblas -lconfig

	/*A
	results = sineRegression_lms(data, f0, fs);
	
	
	printf("\nDados calculados [0]:\n");
	printf("Amplitude: %f\n",gsl_vector_get(results.amplitude,0));
	printf("Fase: %f\n",gsl_vector_get(results.phase_rad,0));
	printf("offset: %f\n\n",gsl_vector_get(results.offset,0));
	
	printf("\nDados calculados [1]:\n");
	printf("Amplitude: %f\n",gsl_vector_get(results.amplitude,1));
	printf("Fase: %f\n",gsl_vector_get(results.phase_rad,1));
	printf("offset: %f\n",gsl_vector_get(results.offset,1));
	* A*/
