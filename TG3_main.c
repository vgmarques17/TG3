#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>
#include "../pmd.h"
#include "../usb-1608FS-Plus.h"
#include "Library/mylib.h"

#define MAX_COUNT     (0xffff)
#define FALSE 0
#define TRUE 1
#define NCHAN	8

int main(int argc, char **argv){
	
	/*Declaraçao de Variáveis*/

	// Variáveis do sistema
	
	libusb_device_handle *udev = NULL; //ponteiro para o dispositivo
	int ret; //parâmetro para saída de funçoes que avaliam o bom funcionamento da placa
	float table_AIN[NGAINS_USB1608FS_PLUS][NCHAN_USB1608FS_PLUS][2]; //Tabela de calibraçao (ranges x canais x 2[inclinaçao e offset ])
	
	//Outras variáveis
	float f_source, f_sampling,time_acq; //Frequência da fonte, frquência de amostragem, tempo de aquisição
	char check; //variável para controle de loops
	int n_ciclos,samples_cicle, N_samples,n_loops,range; /* # total de ciclos por direção, 
												  *	# amostras por ciclo, # amostras total,
												  * número de aquisições, intervalo de medição dos canais */
	int i,m,col = 0; //controle de loops
	uint8_t ranges[8],MUX_door = TRUE, options; //armazena todos os intervalos de medição, porta onde está o MUX, opcoes de aquisicao
	uint8_t channels = 0xff; // Canais a serem usados
	
	FILE *voltage_data;
	char fname_data[100];

	

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
	

	/*Seleçao dos parâmetros da aquisiçao USAR ARQUIVO DE CONFIGURACAO*/ 
	check = 'N';
	do{
		printf("Selecione a frequência da sua fonte (Hz):\n");
		scanf("%f",&f_source);
		printf("Selecione a frequência de amostragem (Hz):\n");
		scanf("%f",&f_sampling);
		/*printf("Por quantos segundos deseja realizar as medidas?\n");
		scanf("%f",&time_acq);*/
		printf("Qual o número total de ciclos da fonte que devem ser aplicados em cada direção?*\n");
		scanf("%d",&n_ciclos);
		
		printf(" Frequência da fonte: %.2f Hz\n Frequência de amostragem: %.2f Hz\n Ciclos em cada direção: %d\n Prosseguir?(Y/N)\n",
				f_source,f_sampling,n_ciclos);
		getchar();
		if((check = getchar())=='Y'||check=='y'){break;};
		
	}while(1);

	printf("\nFazer medidas em qual intervalo de tensão?\n 0: +-10V\n 1: +- 5V\n 3: +-2V\n 5: +-1V\n");
	scanf("%d",&range);
	for(i=0;i<8;i++){
		ranges[i] = range;
	}

	
	

	/*Preparaçao da aquisiçao*/
	samples_cicle = round(f_sampling/f_source); //amostras por ciclo da fonte
	N_samples =  2*n_ciclos*samples_cicle; //número total de amostras
	time_acq = N_samples/f_sampling;
	n_loops = N_samples/samples_cicle; //Número de repetições das medições
	printf("%d amostras serão coletadas para cada canal\nTempo de amostragem: %.2f s\n",N_samples,time_acq);
	
	usbAInScanStop_USB1608FS_Plus(udev); //Para qualquer scan que esteja ocorrendo
	usbAInScanClearFIFO_USB1608FS_Plus(udev);// Limpa o endpoint do FIFO (ie. limpa a fila?)
	usbAInScanConfig_USB1608FS_Plus(udev, ranges);// Configura os ranges de cada canal
	sleep(1);
	
	uint16_t sdataIn[NCHAN*samples_cicle]; //Reserva espaço para os dados de entrada
	float sdataOut[NCHAN*N_samples]; // Reserva espaço para os dados de saída
	
	//Configura o modo de transferência dos dados para o pc
	if (f_sampling < 100.) {
	  options = (IMMEDIATE_TRANSFER_MODE | INTERNAL_PACER_ON);
	} else {
	  options = (BLOCK_TRANSFER_MODE | INTERNAL_PACER_ON);
	}
	
	/*Aquisiçao*/
	
	for(m=0;m<n_loops;m++){
		//1) Manda informaçao para a fonte com MUX
		usbDLatchW_USB1608FS_Plus(udev, MUX_door);
		//2) Coleta dados pelo tempo do frame
		usbAInScanStart_USB1608FS_Plus(udev, samples_cicle, f_sampling, channels, options);
		ret = usbAInScanRead_USB1608FS_Plus(udev, samples_cicle, NCHAN, sdataIn, options);
		
		/*Obs: näo coloquei pra parar pq no teste n tinha, mas pode ser necessário*/
		
		//3) Salva na matriz depois de converter pra volts (Ideal seria fazer isso só no fim) ou multithread!
		for(i=0;i<samples_cicle*NCHAN;i++){
			sdataOut[m*samples_cicle*NCHAN+i] = volts_USB1608FS_Plus(rint(sdataIn[i]*table_AIN[range][(uint8_t)col][0] 
																			 + table_AIN[range][(uint8_t)col][1])
											,range); /*Ajusta valores de acordo com calibracao e converte para Volts*/
			col++;
			if(col==NCHAN){col = 0;}
		}
		//4) Muda direção
		MUX_door = !MUX_door;		
	}
	
/*Como armazenar os dados? Cada vez em uma tabela ou vai aumentando na mesma tabela?*/
	
/*Espaço para demodulaçao*/

	/*Mostrar tabela (temporário)*/
	col=0;
	for(int row=0;row<NCHAN*samples_cicle;row++){
		printf(" %f |",sdataOut[row]);
		col++;	
		if(col==NCHAN){
			printf("\n");
			col=0;
		}
		
	}

	/*Salvar em arquivo .txt (usar funçao)*/
	printf("Deseja salvar? Y/N\n");
	getchar();
	if((check = getchar())=='Y'||check=='y'){
		printf("Escolha um nome para o arquivo .txt:\n");
		scanf("%123s",fname_data);
		strcat(fname_data,".txt");
		voltage_data = fopen(fname_data,"w");
		
		for(i=0;i<NCHAN;i++){
			fprintf(voltage_data,"Channel_%d; ",i );		
		}
		for(int row=0;row<N_samples*NCHAN;row++){
			fprintf(voltage_data," %lf ;",sdataOut[row]);
		}
		col++;
		if(col==NCHAN){col = 0;}
	}

	return 0;
}


/*gcc -g -Wall -I. -o TG3_main TG3_main.c -L. -lmccusb  -lm -L/usr/local/lib -lhidapi-libusb -lusb-1.0 Library/mylib.o*/
