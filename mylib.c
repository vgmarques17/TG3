#include "mylib.h"


/*Intern Functions*/

/*Extern Functions*/

void showTable_daq(float* table, int Nrows, int Ncol){
	int col=0;
	for(int element=0;element<Nrows*Ncol;element++){
		printf(" %f |",table[element]);
		col++;	
		if(col==Ncol){
			printf("\n");
			col=0;
		}
	}
}

void showTable(float** table,int Nrows, int Ncol){
	for(int i=0;i<Nrows;i++){
		for(int j=0;j<Ncol;j++){
			printf("%f |", *((float *)table + (i * Ncol) + j));
		}
		printf("\n");
	}
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

void saveData(float* Data, int Nrows, int Ncol){
	char check;
	FILE *data_file;
	char fname_data[100];
	int col,i,element;
	
	printf("Deseja salvar? Y/N\n");
	getchar();
	
	if((check = getchar())=='Y'||check=='y'){
		printf("Escolha um nome para o arquivo .txt:\n");
		scanf("%123s",fname_data);
		strcat(fname_data,".txt");
		data_file = fopen(fname_data,"w");
		
		for(i=0;i<Ncol;i++){
			fprintf(data_file,"Channel_%d; ",i );		
		}
		for(element=0;element<Nrows*Ncol;element++){
			fprintf(data_file," %lf ;",Data[element]);
		}
		col++;
		if(col==Ncol){col = 0;}
		fclose(data_file);
	}
	
}

int pointer2gsl(gsl_matrix * gsl_M, float* M, int m, int n){
	double element;
	for(int i = 0;i<m;i++){
		for(int j = 0;j<n;j++){
			element = (double)M[i*n+j];
			gsl_matrix_set(gsl_M,i,j,element); 
		}
	}
		
	return 1;
	}
	
int gsl2pointer(float* M, gsl_matrix* gsl_M){
	double element;
	int m = gsl_M->size1, n = gsl_M->size2;
	
	for(int i = 0;i<m;i++){
		for(int j = 0;j<n;j++){
			element = gsl_matrix_get(gsl_M,i,j);
			M[i*n+j] = (float) element; //por causa da precisao da placa, idealmente seria tudo double, mas nao precisa nesse caso
		}
	}
		
	return 1;
	}
	
int gsl2pointer_vector(float* M, gsl_vector* gsl_M){
	double element;
	int m = gsl_M->size;

	for(int i = 0;i<m;i++){
		element = gsl_vector_get(gsl_M,i);
		M[i] = (float) element; //por causa da precisao da placa, idealmente seria tudo double, mas nao precisa nesse caso
	}
	
		
	return 1;
}


int transposeMatrix(float **usr_matrix, int Nrow,int Ncol, float **T_matrix){
	
	for(int row = 0;row<Nrow;row++){
		for(int col = 0;col<Ncol;col++){
			*((float *)T_matrix + (col * Nrow) + row) = *((float *)usr_matrix + (row * Ncol) + col);
		}
	}
	
	return 1;
	
}

int productMatrix(float** A, int rA, int cA,float ** B,int rB, int cB,float** AB){
	/*Allocate memory before using here*/
	
	if(cA!=rB){
		printf("Number of columns in A is different from number of rows in B!\n");
		return -1;
	}
	for(int i=0;i<rA;i++){ //percorre linhas de A e AB
		for(int j=0;j<cB;j++){//percorre colunas de AB
			for(int k=0;k<rB;k++){//percorre linhas de B
				*((float*)AB + (i*cB) + j) += *((float*)A + (i*cA) + k) * *((float*)B + (k*cB) + j) ; 
			}
		}
	}
	return 1;
}

int gsl_productMatrix(gsl_matrix *A, gsl_matrix *B, gsl_matrix * AB){
	double x;
	
	if(A->size2!=B->size1){
		printf("Number of columns in A is different from number of rows in B!\n");
		return -1;
	}
	for(int i=0;i<A->size1;i++){ //percorre linhas de A e AB
		for(int j=0;j<B->size2;j++){//percorre colunas de AB
			for(int k=0;k<B->size1;k++){//percorre linhas de B
				x = gsl_matrix_get(AB,i,j)+gsl_matrix_get(A,i,k)*gsl_matrix_get(B,k,j);
				gsl_matrix_set(AB,i,j,x);
			}
		}
	}
	return 1;
	}


int pinv(gsl_matrix *usr_matrix,gsl_matrix *P_inv){
	
	/*This function calculates the pseudo inverse matrix based on the Moore-Penrose method
	 * using SVD;
	 * The number of rows of the matrix must be >= the number of columns, so that there are*/
	int Nrow = usr_matrix->size1, Ncol = usr_matrix->size2;
	gsl_matrix *U,*Ut,*V,*Sigma_plus,*Sigma,*VS; // SVD's U, U transposed and V matrix, pinv from SVD's Sigma
	gsl_vector *Sigma_array,*work;
	float _temp_value;
		
	/*Initialize necessary Matrices*/
	//U
	U = gsl_matrix_alloc(Nrow,Ncol);
	gsl_matrix_memcpy(U,usr_matrix);
	//Ut
	Ut = gsl_matrix_alloc(Ncol,Nrow);
	//Sigma_array e Sigma
	Sigma_array = gsl_vector_alloc(Ncol);
	Sigma = gsl_matrix_calloc(Ncol,Ncol);
	//Sigma_plus
	Sigma_plus = gsl_matrix_calloc(Ncol,Ncol);
	//V e Vt
	V = gsl_matrix_alloc(Ncol,Ncol);
	//Workspace
	work = gsl_vector_alloc(Ncol);
	VS = gsl_matrix_alloc(Ncol,Ncol);

	/*thin SVD*/
	gsl_linalg_SV_decomp(U,V,Sigma_array,work);
	/* U = m x n
	 * E = n x n
	 * V = n x n*/

	
	/*make Sigma^+*/
	for(int i = 0;i<Ncol;i++){
		_temp_value = gsl_vector_get(Sigma_array,i);
		if(_temp_value<=QUASI_ZERO){
			break;
		}else{
			gsl_matrix_set(Sigma_plus,i,i,1/_temp_value);
		}
		gsl_matrix_set(Sigma,i,i,_temp_value);
	}
	
	/*Transpose U*/
	gsl_matrix_transpose_memcpy(Ut,U);	
	/* get Pseudo inverse
	 * A_plus = V(n*n) * Sigma_plus(n*n) * Ut(n*m)*/
	
	gsl_productMatrix(V,Sigma_plus,VS);
	gsl_productMatrix(VS,Ut,P_inv);
	//gsl2pointer(P_inv,A);
	return 1;
}


fitSine sineRegression_lms(gsl_matrix* Data, float f0,float samp_freq){
	/*Use GSL for all*/
	
	float Ts;
	gsl_vector *time, *omegat; //time vector, angle vector(?)
	gsl_vector *Sines,*Cosines;
	gsl_vector *alpha,*beta,*C;
	gsl_matrix *mat_j, *vec_p, *mat_jI;
	
	gsl_vector *amplitude,*phase_rad;
	gsl_vector *phi_alpha,*phi_beta;
	
	gsl_matrix *gsl_Y;
	float *Y; 
	
	float temp_val_1,temp_val_2;
	gsl_vector *temp_vector_1,* temp_vector_2;
	
	int Nchan,N; //#elements
	fitSine results;
	
	//Initialization
	
	Ts = 1/samp_freq;
	Nchan = Data->size2;
	N = Data->size1; //size of data		
	
	time = gsl_vector_calloc(N); //alloc memory for arrays
	omegat = gsl_vector_calloc(N);
	
	mat_j = gsl_matrix_alloc(N,3); 
	mat_jI = gsl_matrix_alloc(3,N);
	
	vec_p = gsl_matrix_alloc(3,Nchan); 
	
	alpha = gsl_vector_alloc(Nchan);
	beta = gsl_vector_alloc(Nchan);
	C =  gsl_vector_alloc(Nchan);
	
	amplitude = gsl_vector_alloc(Nchan);
	phase_rad = gsl_vector_alloc(Nchan);
	phi_alpha = gsl_vector_alloc(Nchan);
	phi_beta = gsl_vector_alloc(Nchan);
	
	gsl_Y = gsl_matrix_alloc(N,Nchan);
	Y = malloc(N*Nchan*sizeof(float));
	
	//final results
	results.amplitude = malloc((Nchan*sizeof(float)));
	results.phase_rad = malloc((Nchan*sizeof(float)));
	results.offset = malloc((Nchan*sizeof(float)));
	results.alpha = malloc((Nchan*sizeof(float)));
	results.beta = malloc((Nchan*sizeof(float)));
	results.omegat = malloc((Nchan*sizeof(float)));
	results.y = malloc((N*Nchan*sizeof(float)));
	
	
	//temporary
	temp_vector_1 = gsl_vector_alloc(N);
	temp_vector_2 = gsl_vector_alloc(N); //livrar vetores dps em todas as fun¢øes
	
	//--------------------//
	for(int i = 0; i<N;i++){gsl_vector_set(time,i,(double)i*Ts);} //creates time vector
	
	for(int i = 0; i<N;i++){gsl_vector_set(omegat,i,2*PI*f0*gsl_vector_get(time,i));}//create omega vector;
	
	/*Compute mat_J*/
	for(int i=0;i<N;i++){
		temp_val_1 = gsl_vector_get(omegat,i);
		gsl_matrix_set(mat_j,i,0,(float)sin(temp_val_1)); //first column  = sines
		gsl_matrix_set(mat_j,i,1,(float)cos(temp_val_1)); // second columns = cosines
		gsl_matrix_set(mat_j,i,2,1.0);					 // third column = 1
	}

	pinv(mat_j,mat_jI);
	
	gsl_productMatrix(mat_jI,Data,vec_p);//nessa matriz estao contidos alfa, beta e C para cada canal que houver em Data

	//separa vec_p
	gsl_matrix_get_row(alpha,vec_p, 0);
	gsl_matrix_get_row(beta,vec_p, 1);
	gsl_matrix_get_row(C,vec_p, 2);

	//amplitude e fase
	for(int i = 0;i<Nchan;i++){
		temp_val_1 = sqrt(pow(gsl_vector_get(alpha,i),2)+pow(gsl_vector_get(beta,i),2));
		gsl_vector_set(amplitude,i,temp_val_1);
		
		temp_val_2 = atan2(gsl_vector_get(beta,i),gsl_vector_get(alpha,i));
		gsl_vector_set(phase_rad,i,temp_val_2);
				
	}
	//Phi alfa e phi beta
	phi_alpha = alpha;
	gsl_vector_div(phi_alpha,amplitude);
	
	phi_alpha = beta;
	gsl_vector_div(phi_beta,amplitude);
	
	for(int i = 0;i<Nchan;i++){
		temp_val_1 = acos(gsl_vector_get(phi_alpha,i));
		gsl_vector_set(phi_alpha,i,temp_val_1);
		
		temp_val_2 = acos(gsl_vector_get(phi_beta,i));
		gsl_vector_set(phi_beta,i,temp_val_2);	
	}
	
	//Y
	for(int i = 0;i<N;i++){
		temp_vector_1 = alpha;
		temp_vector_2 = beta;
		gsl_vector_scale(temp_vector_1,sin(2*PI*f0*gsl_vector_get(time,i)));
		gsl_vector_scale(temp_vector_2,cos(2*PI*f0*gsl_vector_get(time,i)));
		
		gsl_vector_add(temp_vector_1,temp_vector_2);
		gsl_vector_add(temp_vector_1,C);
		gsl_matrix_set_row(gsl_Y,i,temp_vector_1);				
	}
	
	// Mandar para os resultados, rever necessidade de usar ponteiros no programa principal
	gsl2pointer_vector(results.amplitude,amplitude);
	gsl2pointer_vector(results.phase_rad,phase_rad);
	gsl2pointer_vector(results.offset,C);
	gsl2pointer_vector(results.alpha,alpha);
	gsl2pointer_vector(results.beta,beta);
	gsl2pointer_vector(results.omegat,omegat);
	gsl2pointer(results.y,gsl_Y);
	
	
	return results;
}
/* ycos(k,:) = beta.*cos(2*pi*f0*discrete_time(k));
   ysin(k,:) = alpha.*sin(2*pi*f0*discrete_time(k));
   y(k,:) = ysin(k,:) + ycos(k,:) + C; 
 * */


/* Nota para o futuro: se quiser passar uma matriz como parâmetro, 
 * utilizar as dimensoes como parâmetros ANTES dela, ou mandar como
 * ponteiro
 * 
 * APRENDER A USAR PONTEIROS DIREITO, ESSA BIBLIOTECA FUNCIONA SÓ COM TAMANHOS PRÉ DEFINIDOS DE MATRIZ*/
