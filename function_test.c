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
	
	int samples_cicle,n_cicles;
	double f0,fs,samples,time;


	gsl_matrix * data, *electrode_pairing;

	config_t cfg_params;
	config_setting_t *setting;
	char cfg_file[100],txt_file[100];
	
	Sys_Results Output;
	
	
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
	
	config_lookup_float(&cfg_params,"f_sampling",&fs);
	config_lookup_float(&cfg_params,"f_source",&f0);
	config_lookup_float(&cfg_params,"acq_time",&time);
	config_lookup_int(&cfg_params,"times_direction",&n_cicles);
	samples = fs*time;
	samples_cicle = round(samples/n_cicles);

	
	//matriz com pareamento dos eletrodos
	// Duas primeiras linhas: DDP músculo (L,T)
	// Duas últimas linhas: DDP sentinela (L,T)
	
	electrode_pairing = gsl_matrix_alloc(4,2);
	
	gsl_matrix_set(electrode_pairing,0,0,1); gsl_matrix_set(electrode_pairing,0,1,0); // [1 0] 
	gsl_matrix_set(electrode_pairing,1,0,2); gsl_matrix_set(electrode_pairing,1,1,3); // [2 3]
	gsl_matrix_set(electrode_pairing,2,0,4); gsl_matrix_set(electrode_pairing,2,1,5); // [4 5]
	gsl_matrix_set(electrode_pairing,3,0,6); gsl_matrix_set(electrode_pairing,3,1,7); // [6 7]
	
	//
	printf("Enter .txt file name/path:\n");
	scanf("%123s",txt_file);
	strcat(txt_file,".txt");
	
	data = readFile_gsl(txt_file,(int) samples,8,TRUE);
	

	Output = calculateImpedance(data, samples_cicle, electrode_pairing, f0, fs);
	saveFile_gsl(cfg_params,Output.impedance_data, 1);
	saveFile_gsl(cfg_params,Output.phasor_data, 2);

	
	
	for(int i = 0;i<Output.phasor_data->size1;i++){
		for(int j = 0;j<Output.phasor_data->size2;j++){
			printf("%0.4f| ",gsl_matrix_get(Output.phasor_data,i,j));
		}
		printf("\n");
	}
	printf("\n\n");
	for(int i = 0;i<Output.impedance_data->size1;i++){
		for(int j = 0;j<Output.impedance_data->size2;j++){
			printf("%0.4f| ",gsl_matrix_get(Output.impedance_data,i,j));
		}
		printf("\n");
	}
	
	
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
	
	/*B	
	
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
	
	impedance = gsl_matrix_alloc(factor,12);*/
