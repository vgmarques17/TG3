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
	
	/*Declaraçao de Variáveis*/

	// Variáveis do sistema
	
	libusb_device_handle *udev = NULL; //ponteiro para o dispositivo
	int ret; //parâmetro para saída de funçoes que avaliam o bom funcionamento da placa
	float table_AIN[NGAINS_USB1608FS_PLUS][NCHAN_USB1608FS_PLUS][2]; //Tabela de calibraçao (ranges x canais x 2[inclinaçao e offset ])
	
	//Outras variáveis
	config_t cfg_params;
	char cfg_file[100];
	
	double f_source, f_sampling,acq_time,T_s; //Frequência da fonte, frquência de amostragem, tempo de aquisição, periodo de aquisicao
	
	int samples_cicle, N_samples,n_loops,range, cicles_time,times_direction,N_chan; /* # amostras por ciclo, # amostras total,
														* número de aquisições, #ciclos por dire¢ao, #numero de aplicacoes por direcao,
														* Numero de canais */
	int col = 0; //controle de loops
	uint8_t ranges[8],MUX_door = TRUE, options; //armazena todos os intervalos de medição, porta onde está o MUX, opcoes de aquisicao
	uint8_t channels = 0xff; // Canais a serem usados
	
	gsl_matrix *Data_Out,*electrode_pairing;
	
	Sys_Results Output;

	

	/*Inicializaçao da placa*/
	udev = NULL;
	ret = libusb_init(NULL);
	if(ret<0){
		perror("Falha ao iniciar a libusb!");
		exit(1);
	}	
	if((udev = usb_device_find_USB_MCC(USB1608FS_PLUS_PID,NULL))){
		printf("Dispositivo encontrado!\n");
	}else{
		printf("Falha, placa USB 1608FS-Plus não encontrada!\n");
		return 0;
	}
	
	usbBuildGainTable_USB1608FS_Plus(udev, table_AIN); //obtém parâmetros de calibração
	
	usbDTristateW_USB1608FS_Plus(udev,0xfe); /*Configura quais portas digitais são
											  *entradas ou saídas. 0 = saída, 1 = entrada
											  *Ler como binário; Configuração do MUX*/
											  
	printf("\n-----------------MIOGRAFIA POR IMPEDÂNCIA ELÉTRICA----------------- \n\n ");
	

	/*Parâmetros da aquisiçao USAR ARQUIVO DE CONFIGURACAO*/ 
	
	
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
	
	/*Decision Tree para determinar tipo de configura¢ao*/

	if(config_lookup_float(&cfg_params,"f_sampling",&f_sampling)==0.0){
		//Configura¢ao usando periodo
		config_lookup_float(&cfg_params,"T_s",&T_s);
		f_sampling = 1/T_s;

	}else{
		//Configura¢ao usando frequencia de amostragem
		T_s = 1/f_sampling;
	}
	
	config_lookup_float(&cfg_params,"f_source",&f_source); // frequencia da fonte
		
	if(config_lookup_float(&cfg_params,"acq_time",&acq_time) != 0){
		
		N_samples = acq_time*f_sampling;
		config_lookup_int(&cfg_params,"cicles_time",&cicles_time); //Da erro se colocado direto no if
		
		if(cicles_time == 0){
			//Configura¢æo com times_direction e acq_time
			config_lookup_int(&cfg_params,"times_direction",&times_direction);
			samples_cicle = (int)round(N_samples/times_direction);				
			cicles_time = round(N_samples/(samples_cicle*2));
		}else{
			//Configura¢ao com cicles_time e acq_time
			samples_cicle = round(N_samples/(cicles_time*2));
			times_direction = round(N_samples/samples_cicle);
		}
	}else{
		//Configura¢ao com cicles_time e times_direction
		config_lookup_int(&cfg_params,"cicles_time",&cicles_time);
		config_lookup_int(&cfg_params,"times_direction",&times_direction);
		N_samples = times_direction*2*(cicles_time*f_sampling/f_source);
		samples_cicle = round(N_samples/(cicles_time*2));
		acq_time = N_samples*T_s;
	}
	/*colocar verifica¢øes para erros*/
	

	
	printf("\n\nSumário da Aquisiçao \n----------------------\n\n");
	
	printf("Frequência de Amostragem [Hz]: %0.2f \n",f_sampling);
	printf("Período de Amostragem [ms]: %0.2f \n",T_s*1000);
	printf("Frequência da Fonte [Hz]: %0.2f \n",f_source);
	printf("Número total de amostras: %i \n",N_samples);
	printf("Ciclos da fonte para cada aplicacao em cada direcao: %i \n",cicles_time);
	printf("Números de aplicacoes em cada direcao: %i \n",times_direction);
	printf("Número de amostras coletadas em cada aplicacao: %i \n",samples_cicle);
	printf("Tempo estimado de aquisicao [s]: %0.2f \n",acq_time);
	
	//Colocar confirmaçao do usuario
	
	
	/*Preparaçao da aquisiçao*/
	
	n_loops = N_samples/samples_cicle; //Número de repetições das medições
	config_lookup_int(&cfg_params,"range",&range);
	config_lookup_int(&cfg_params,"N_chan",&N_chan);
	for(int i =0;i<N_chan;i++){ranges[i] = range;}
	
	usbAInScanStop_USB1608FS_Plus(udev); //Para qualquer scan que esteja ocorrendo
	usbAInScanClearFIFO_USB1608FS_Plus(udev);// Limpa o endpoint do FIFO (ie. limpa a fila?)
	usbAInScanConfig_USB1608FS_Plus(udev, ranges);// Configura os ranges de cada canal
	sleep(1);

	uint16_t sdataIn[N_chan*samples_cicle]; //Reserva espaço para os dados de entrada
	float sdataOut[N_chan*N_samples]; // Reserva espaço para os dados de saída
	Data_Out = gsl_matrix_alloc(N_samples,N_chan); // Reserva espaço para os dados de saída em GSL
	
	
	//Configura o modo de transferência dos dados para o pc
	if (f_sampling < 100.) {
	  options = (IMMEDIATE_TRANSFER_MODE | INTERNAL_PACER_ON); //problema nessa opcao
	} else {
	  options = (BLOCK_TRANSFER_MODE | INTERNAL_PACER_ON);
	}
	
	
	/*Aquisiçao*/
	printf("\nInício da Aquisiçao\n");
	for(int m=0;m<n_loops;m++){
		//1) Manda informaçao para a fonte com MUX
		usbDLatchW_USB1608FS_Plus(udev, MUX_door);
		//2) Coleta dados pelo tempo do frame
		usbAInScanStart_USB1608FS_Plus(udev, samples_cicle, f_sampling, channels, options);
		ret = usbAInScanRead_USB1608FS_Plus(udev, samples_cicle, N_chan, sdataIn, options);
		
		/*Obs: näo coloquei pra parar pq no teste n tinha, mas pode ser necessário*/
		
		//3) Salva na matriz depois de converter pra volts (Ideal seria fazer isso só no fim) ou multithread! 
				// Outra opcao seria só copiar a tabela para a tabela maior,e depois converter. Deve gastar menos tempo
		for(int i = 0;i<samples_cicle*N_chan;i++){
			sdataOut[m*samples_cicle*N_chan+i] = volts_USB1608FS_Plus(
													rint(sdataIn[i]*table_AIN[range][(uint8_t)col][0] 
													+ table_AIN[range][(uint8_t)col][1]),range
												); /*Ajusta valores de acordo com calibracao e converte para Volts*/
			col++;
			if(col==N_chan){col = 0;}
		}
		//4) Muda direção
		MUX_door = !MUX_door;		
	}
	printf("\nFim da Aquisiçao\n");
	
	//showTable_1Df(sdataOut,100,8);
	
	//matriz com pareamento dos eletrodos
	// Duas primeiras linhas: DDP músculo (L,T)
	// Duas últimas linhas: DDP sentinela (L,T)
	
	electrode_pairing = gsl_matrix_alloc(4,2);
	
	gsl_matrix_set(electrode_pairing,0,0,2); gsl_matrix_set(electrode_pairing,0,1,1); //  VL
	gsl_matrix_set(electrode_pairing,1,0,0); gsl_matrix_set(electrode_pairing,1,1,3); //  VT
	gsl_matrix_set(electrode_pairing,2,0,6); gsl_matrix_set(electrode_pairing,2,1,7); //  IL
	gsl_matrix_set(electrode_pairing,3,0,4); gsl_matrix_set(electrode_pairing,3,1,5); //  IT
	
	
	/*Espaço para demodulaçao*/
	pointer2gsl_matrix(Data_Out,sdataOut,N_samples,N_chan);
	
	Output = calculateImpedance(Data_Out,samples_cicle,electrode_pairing,f_source,f_sampling);

	/*Salvar em arquivo .txt (usar funçao)*/
	//ver se usuario quer salvar
	
	saveFile_gsl(cfg_params,Data_Out, 0);	
	saveFile_gsl(cfg_params,Output.impedance_data, 1);
	saveFile_gsl(cfg_params,Output.phasor_data, 2);

	
	config_destroy (&cfg_params);	

	return 0;
}


/*gcc -g -Wall -I. -o TG3_main TG3_main.c -L. Library/mylib.o -lmccusb -L/usr/local/lib -lhidapi-libusb -lusb-1.0 -lm -lgsl -lgslcblas -lconfig */
