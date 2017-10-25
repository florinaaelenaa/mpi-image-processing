/*
 * Barbu Florina
 * Procesarea de imagini folosind rețele MPI
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi/mpi.h"


/**** TAGURI *******************************************************************************************************************************************************/
#define PING 	   1
#define ECHO 	   2
#define TO_PROCESS 3
#define PROCESSED  4
#define DIM 	   5
#define STATISCS   6
#define DONE	   7
/******************************************************************************************************************************************************************/

int** aloc_matrix(int n, int m);

int** matrix_border(int rows, int cols, int **mat);

int** read_photo(int *n, int *m, FILE *f, char *pgm_type, int *max_px_val);

void display_photo(int n, int m, char* pgm_type, int max_px_val, int **matrix, FILE *f);

int** create_filter(char *filter_name, int *factor);

int** apply_filter(int rank, int** filter, int no_processes, int *lines_processed, int rows, int cols, int ** picture, int factor);

void create_new_pic(int rank, char *file_name, int no_processes, char *pgm_type, int max_read);

int main(int argc, char** argv) 
{

/****** DATE ***********************************************************************************************************************************************************/
	int **graph;
	int **copy;
	int **void_t;
	int matrix_buff[100][100] = {0}; //buffer in care citesc topologia inainte s-o pun in graph pe care il aloc dinamic pt a fi la fix cat imi trebuie
	char buffer[101]; //buffer in care citesc cate o linie din fisier pt a-mi crea topologia

	int echoCount = 0; // cate ecouri trebuie sa primesc
	int children = 0;  // nr de copii ai ficarui nod
	int dim_to_send = 0; //nr de linii pe care fiecare proces trebuie sa-l trimita mai departe copiilor
	volatile int dim_received = 0; //nr de linii pe care le primeste fiecare proces fie pt procesare/ fie sa le distribuie mai departe
	int dim;
	int *dim_vect; 
	int lines_processed = 0;
	int proba = 1;
	int tag;
	int parent = 0;
	int source = 0;
	int num_nodes = 0;
	char * pch;
	
	int i, j, citit = 0, index;
	int no_pics;																	// poze de procesat in fisierul imagini.in
	int n; 																			// nr de linii
	int m; 																			// nr de coloane
	int **picture;
	int rank, nProcesses;

	char pgm_type[3];																/* No  Colors				   File-extesion
																					   P1: 0–1   (black & white)   .pbm
																					   P2: 0–255 (gray scale)      .pgm
																					   P3: 0–255 (RGB)        	   .ppm
																					*/	
	int max_px_val;																	//valoarea maxima care poate fi luata de un px din poza

	char *p_file, *t_file, *s_file;
	char p_name[30], o_name[30];													//numele imaginii ce urmeaza a fi procesata si numele imaginii dupa procesare												
	char filter_name[15];															//tipul de filtru citit de pe fiecare linie din imagini.in
	int **filter;
	int factor;
	
	FILE *pics_file, *topology_file, *statistics_file, *pic_name, *out_file;
/********************************** *******************************************************************************************************************************/


	// Initializare mediu MPI
    n = MPI_Init(&argc, &argv);
    MPI_Status status;
    
    if (n != MPI_SUCCESS) 
    {
    	printf ("Error starting MPI program. Terminating.\n");
    	MPI_Abort(MPI_COMM_WORLD, n);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
   
	t_file = strdup(argv[1]);													//imagini.in
	p_file = strdup(argv[2]);													//topologie.in
	s_file = strdup(argv[3]);													//statistica.out

	pics_file       = fopen(p_file, "rt");										//pointer la fisierul imagini.in
	topology_file   = fopen(t_file, "rt");										//pointer la fisierul topologie.in
	statistics_file = fopen(s_file, "wt");										//pointer la fisierul statistica.out

	fscanf(pics_file, "%d", &no_pics);											//acest numar de poze va trebui decrementat odata ce s-a terminat procesarea unei poze

	

/*** Imi construiesc topologia ca matrice de adiacenta ************************************************************************************************************/
	while( !feof(topology_file) )
	{
		
		fgets(buffer, 100, topology_file);
		
		pch = strtok (buffer," ");
		
		int count = 0; //ma ajuta sa sar peste primul string in care am nodul curent, ai carui vecini urmeaza

		while (pch != NULL)
		{
			index = atoi(pch);
			if(count > 0)
				matrix_buff[num_nodes][index] = 1;
			
			pch = strtok (NULL, " ");
			
			count ++;		
		}
		num_nodes++;	
	}

	graph = aloc_matrix(num_nodes, num_nodes);									//topologia curenta a procesului
	copy  = aloc_matrix(num_nodes, num_nodes);
	void_t = aloc_matrix(num_nodes, num_nodes); 								//o folosim ca sa eliminam ciclurile
/********************************** *******************************************************************************************************************************/



/*** Calculeaza arbore de acoperire *******************************************************************************************************************************/
	for(i = 0; i < num_nodes; ++i)
	{	
		for (j = 0; j < num_nodes; ++j)
		{
			graph[i][j] = matrix_buff[i][j];
		}
		memset(copy[i],   0, num_nodes * sizeof(int));
		memset(void_t[i], 0, num_nodes * sizeof(int));
	}
	
	//vrem ca un nod sa-si cunoasca doar vecinii
	for(i = 0; i < num_nodes; i++)
		for(j = 0; j < num_nodes; j++)
			if(rank != i)
				graph[i][j] = 0; //deci facem 0 in rest in matricea de adiacenta


	MPI_Barrier(MPI_COMM_WORLD);
	if(rank != 0)
	{
		
		MPI_Recv(&(copy[0][0]), num_nodes * num_nodes, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		parent = status.MPI_SOURCE; //aflu parintii proceselor care sunt diferite de cel initiator

	}
	for(i = 0; i < nProcesses; i++) 
	{
		if(graph[rank][i] && i != parent) 
		{
			MPI_Send(&proba, 1, MPI_INT, i, PING, MPI_COMM_WORLD); //send dummy message => ping
			echoCount++;										   //la finalul for-ului aflu practic cati copii are fiecare nod
		}
		children = echoCount;
	}
	

	while(echoCount > 0) 
	{
		
		MPI_Recv(&(copy[0][0]), num_nodes * num_nodes, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);// ping sau echo

		tag = status.MPI_TAG;
		source = status.MPI_SOURCE;
		
		if(tag == PING)
		{
			MPI_Send(&(void_t[0][0]), num_nodes * num_nodes, MPI_INT, source, ECHO, MPI_COMM_WORLD);	//daca exista cicluri trimitem topologie vida ca sa stim sa stergem legaturile
		}
		else if(tag == ECHO) 
		{
			int schimbat = 0;
			
			for(i = 0; i < num_nodes; i++)
				for(j = 0; j < num_nodes; j++)
				{
					int old = graph[i][j];

					graph[i][j] = graph[i][j] || copy[i][j]; // face Merge cu valorile noi
					
					if(graph[i][j] != old)
						schimbat = 1;
				}

			echoCount--;
			
			if(schimbat == 0) 
			{
				graph[rank][source] = 0; //eliminam ciclurile
				children --;
			}
		}
	}

	if(rank != 0) 
	{
		MPI_Send(&(graph[0][0]), num_nodes * num_nodes, MPI_INT, parent, ECHO, MPI_COMM_WORLD);
		
	}
/******************************************************************************************************************************************************************/

	MPI_Barrier(MPI_COMM_WORLD);												//toate proceselor asteapta ca procesul MASTER sa citeasca din imagini.in si sa face bcast-ul




/*** Incep procesarea pozelor *************************************************************************************************************************************/
	while(no_pics >= 0)
    {
    	if(no_pics > 0)
    	{
    		filter = (int**)malloc(3 * sizeof(int*));
		    int *buf2 = (int*)malloc(9 * sizeof(int));
			for(i = 0; i < 3; i++) 
				filter[i] = &(buf2[i * 3]);

		    if(rank == 0)
		    {
				fscanf(pics_file, "%s", filter_name);								//tipul filtrului ce urmeaza a fi aplicat pe imaginea citita
				fscanf(pics_file, "%s", p_name);									//string: nume_imagine.pgm
		    	fscanf(pics_file, "%s", o_name);
		    	
		    	pic_name  = fopen(p_name, "rt");
		    	out_file  = fopen(o_name, "wt");

		    	filter = create_filter(filter_name, &factor);

		    	picture = read_photo(&n, &m, pic_name, pgm_type, &max_px_val);
		    	n += 2;
		    	m += 2;
		    }

		    MPI_Bcast(&(filter[0][0]), 9, MPI_INT, 0, MPI_COMM_WORLD);
		    MPI_Bcast(&factor, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD); 
			
			if(rank == 0)
		    {
		    	//am citit poza si acum vrem s-o distribuim in mod egal proceselor copii
		    	//pentru asta calculam cati copii are radacina pentru a vedea cum impartim poza pt a o trimite la procesat
		    	dim_to_send = n / children;

		    	// TRIMIT LA COPII
		 		for(i = 0, j = 0; i < num_nodes; i++)
		 		{
		 			if(graph[0][i] == 1 && i != parent) //daca exista legatura in matricea de adiacenta e copil => fac send
		 			{	

		 				int fromLine = j * dim_to_send - (j == 0 ? 0 : 1);
		 				int toLine = (j + 1) * dim_to_send;
		 				if(j == children - 1) // last children
		 				{
		 					toLine = (j + 2) * dim_to_send;
		 				}

		 				toLine = toLine >= n ? n - 1 : toLine;
		 				int linesCount = toLine - fromLine + 1;

		 				printf("rank %d> Send to %d from line = %d to %d, %d lines\n", rank, i, fromLine, toLine, linesCount);

		 				MPI_Send(&linesCount, 1, MPI_INT, i, DIM, MPI_COMM_WORLD); //ii trimit intai cate linii sa se astepte sa primeasca
		 				MPI_Send(&(picture[fromLine][0]), linesCount * m, MPI_INT, i, TO_PROCESS, MPI_COMM_WORLD);

		 				j++;
		 			}

		 		}
		 			
		 		// PRIMESC DE LA COPII FIX CAT AM TRIMIS
		 		for(i = 0, j = 0; i < num_nodes; i++)
		 		{
		 			if(graph[rank][i] == 1 && i != parent) //daca exista legatura in matricea de adiacenta e copil => fac send
		 			{
		 				int fromLine = j * dim_to_send - (j == 0 ? 0 : 1);
		 				int toLine = (j + 1) * dim_to_send;

		 				if(j == children - 1) 
		 				{
		 					toLine = (j + 2) * dim_to_send;
		 				}

		 				toLine = toLine >= n ? n - 1 : toLine;
		 				int linesCount = toLine - fromLine + 1;

		 				printf("rank %d> Recv from %d from line = %d to %d, %d lines\n", rank, i, fromLine, toLine, linesCount);
		 				MPI_Recv(&(picture[fromLine+1][0]), (linesCount - 2) * m, MPI_INT, i, PROCESSED, MPI_COMM_WORLD, &status);

		 				j++;

		 			}
				}

		     	n -= 2;
		    	display_photo(n, m, pgm_type, max_px_val, picture, out_file); 
			}
			else
			{
				MPI_Recv(&dim, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //facem recieve cu nr de linii pe care ar trebui sa-l primim
		 		printf("rank %d> recv %d lines from %d\n", rank, dim, status.MPI_SOURCE);
		 		picture = aloc_matrix(dim, m); //ca proces intermediar/frunza imi creez si eu un buffer picture, in care primesc datele de la parinte
		 		
		 		MPI_Recv(&(picture[0][0]), dim * m, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		 	
		 		if(children > 0) //verific daca sunt proces intermediar
				{

/************* Trimit datele catre copii ca proces intermediar **************************************************************************************************/
					dim_to_send = dim / children; //numarul de linii impartit la cati copii are nodul

			 		for(i = 0, j = 0; i < num_nodes; i++)
			 		{
			 			if(graph[rank][i] == 1 && i != parent) //daca exista legatura in matricea de adiacenta e copil => fac send
			 			{
			 				int fromLine = j * dim_to_send - (j == 0 ? 0 : 1);
			 				int toLine = (j + 1) * dim_to_send;
			 				if(j == children - 1) 
			 				{
			 					toLine = (j + 2) * dim_to_send;
			 				}

			 				toLine = toLine >= dim ? dim - 1 : toLine;
			 				int linesCount = toLine - fromLine + 1;

			 				printf("rank %d> Send to %d from line = %d to %d, %d lines\n", rank, i, fromLine, toLine, linesCount);

			 				MPI_Send(&linesCount, 1, MPI_INT, i, DIM, MPI_COMM_WORLD); //ii trimit intai cate linii sa se astepte sa primeasca
			 				MPI_Send(&(picture[fromLine][0]), linesCount * m, MPI_INT, i, TO_PROCESS, MPI_COMM_WORLD);

			 				j++;

			 			}
					}
/******************************************************************************************************************************************************************/
		

/************* Astept datele procesate de copii, ca proces intermediar *********************************************************************************************/
				
					for(i = 0, j = 0; i < num_nodes; i++)
			 		{
			 			if(graph[rank][i] == 1 && i != parent) //daca exista legatura in matricea de adiacenta e copil => fac send
			 			{
			 				int fromLine = j * dim_to_send - (j == 0 ? 0 : 1);
			 				int toLine = (j + 1) * dim_to_send;
			 				if(j == children - 1) 
			 				{
			 					toLine = (j + 2) * dim_to_send;
			 				}

			 				toLine = toLine >= dim ? dim - 1 : toLine;
			 				int linesCount = toLine - fromLine + 1;

			 				printf("rank %d> Recv from %d from line = %d to %d, %d lines\n", rank, i, fromLine, toLine, linesCount);

			 				//MPI_Send(&linesCount, 1, MPI_INT, i, DIM, MPI_COMM_WORLD); //ii trimit intai cate linii sa se astepte sa primeasca
			 				MPI_Recv(&(picture[fromLine + 1][0]), (linesCount-2) * m, MPI_INT, i, PROCESSED, MPI_COMM_WORLD, &status);

			 				j++;

			 			}
					}

/******************************************************************************************************************************************************************/
	


/*************** Trimit poza procesata la parinte ca si proces intermediar ****************************************************************************************/				
					MPI_Send(&(picture[1][0]), (dim - 2) * m, MPI_INT, parent, PROCESSED, MPI_COMM_WORLD); //cu taggul ca datele sunt procesate
					
/******************************************************************************************************************************************************************/

				}

				else if(children == 0) //inseamna ca sunt pe nodurile frunza care fac prelucrarea
				{
/*************** Procesez poza ca si proces copil *****************************************************************************************************************/
					picture = apply_filter(rank, filter, num_nodes, &lines_processed, dim, m, picture, factor);

/******************************************************************************************************************************************************************/


/*************** Trimit poza procesata la parinte ca si proces copil **********************************************************************************************/				
					
					MPI_Send(&(picture[1][0]), (dim - 2) * m, MPI_INT, parent, PROCESSED, MPI_COMM_WORLD); //cu taggul ca datele sunt procesate
/******************************************************************************************************************************************************************/
					
				}
					 		
			}
/******************************************************************************************************************************************************************/
		}
		else if(no_pics == 0)
		{
			if(rank == 0)
			{
/*************** Trimit tag-ul de terminare de procesat poze ca proces radacina ************************************************************************************/		
				
				for(i = 0; i < num_nodes; i++) 
					if(graph[rank][i] && i != parent) 
						MPI_Send(&proba, 1, MPI_INT, i, DONE, MPI_COMM_WORLD);

/******************************************************************************************************************************************************************/

/*************** Astept statitisca sa o scriu in fisier  ***********************************************************************************************************/				
				
				for(i = 0; i < num_nodes; i++) 
					if(graph[rank][i] == 1 && i != parent) 
					{
						int *buf =(int*)calloc(num_nodes, sizeof(int));
						MPI_Recv(&(buf[0]) , num_nodes, MPI_INT, i, STATISCS, MPI_COMM_WORLD, &status);
						
						for(j = 0; j < num_nodes; j++)
							if(buf[j] > 1)
								graph[rank][j] = buf[j];
					}
					
				for(i = 0; i < num_nodes; i++) 
				{	
					if(graph[rank][i] > 1)	
						fprintf(statistics_file, "%d: %d\n", i, graph[rank][i]);											
						
					else 
						fprintf(statistics_file, "%d: 0\n", i);											
				}			

/******************************************************************************************************************************************************************/
			}
			else
			{
/*************** Astept tagul de terminare si trimit mai departe daca sunt proces intermediar, altfel ca frunza trimit statistica *********************************/				
				
				MPI_Recv(&proba, 1, MPI_INT, MPI_ANY_SOURCE, DONE, MPI_COMM_WORLD, &status);
				if(children > 0)
				{	
					for(i = 0; i < num_nodes; i++)  
						if(graph[rank][i]  && i != parent) 
							MPI_Send(&proba, 1, MPI_INT, i, DONE, MPI_COMM_WORLD);  //trimite DONE copiilor


					for(i = 0; i < num_nodes; i++) 
						if(graph[rank][i] == 1 && i != parent) 
						{
							int *buf =(int*)calloc(num_nodes, sizeof(int));

							MPI_Recv(&(buf[0]) , num_nodes, MPI_INT, i, STATISCS, MPI_COMM_WORLD, &status); //primeste statistica de la copii
							
							for(j = 0; j < num_nodes; j++)
								if(buf[j] > 1)	
									graph[rank][j] = buf[j];
							
						}	

					MPI_Send(&(graph[rank][0]), num_nodes, MPI_INT, parent, STATISCS, MPI_COMM_WORLD);	 //trimite statistica la parinti
					
					
				}
				else
				{	
					printf("rank %d> lines_processed: %d\n", rank, lines_processed);
					graph[rank][rank] = lines_processed;
					
					MPI_Send(&(graph[rank][0]), num_nodes, MPI_INT, parent, STATISCS, MPI_COMM_WORLD);
				}

/******************************************************************************************************************************************************************/
			
			}
		}
			
		MPI_Barrier(MPI_COMM_WORLD);
		no_pics--;
		
	}
    
    MPI_Finalize();
	
	return 0;
}

int** aloc_matrix(int rows, int cols)
{
	int **mat;
	int i,j;
	
	mat = (int **)malloc(rows * sizeof(int*));

	int *buf = (int*)malloc(rows * cols * sizeof(int));
	for(i = 0; i < rows; i++) 
		mat[i] = &(buf[i * cols]);

    return mat;
}

int** matrix_border(int rows, int cols, int **mat)
{
	int i,j;

	//Bordam matricea cu 0 ca sa putem aplica ulterior filtrele
	for (i = 0; i < rows; i++)
	{
		for (j = 0; j < cols; j++)
		{
			if( i == 0 || i == (rows - 1) || j == 0 || j == (cols - 1))
				mat[i][j] = 0;
		}
	}
   
    return mat;
}


int** read_photo(int *n, int *m, FILE *f, char *pgm_type, int *max_px_val)
{
	int **matrix;
	char aux[201];
	int i, j;
		
	//Citire tipul imaginii
    fgets(pgm_type, 4, f);
    printf("%s", pgm_type);
    
    //Citim comentariul
    fgets(aux, 200, f);
    puts(aux);

    //Citire dimensiuni imagine
    fscanf(f, "%d", m);
    fscanf(f, "%d", n);
    fscanf(f, "%d", max_px_val);
	
	matrix = aloc_matrix(*n + 2, *m + 2);
	matrix = matrix_border(*n + 2, *m + 2, matrix);

    for(i = 1; i <= *n; i++)
        for(j = 1; j <= *m; j++)
        {
            if(j < *m - 1)
                fscanf(f, "%d ", &matrix[i][j]); 
            else
                fscanf(f, "%d\n", &matrix[i][j]);
            
        }
  

	return matrix;
}

void display_photo(int n, int m, char* pgm_type, int max_px_val, int **matrix, FILE *f)
{
	int i, j;
	char *f_new;
	
	fprintf(f, "%s\n", pgm_type);												//tipul fisierului (pgm)
	fprintf(f, "%d %d\n", m, n);												//dimensiunea (n x m) a pozei
	fprintf(f, "%d\n", max_px_val);												//valoarea maxima a unui px: 255
    
    for(i = 0; i < n; i++)
        for(j = 1; j <= m; j++)
        	fprintf(f, "%d\n", matrix[i][j]);									//fiecare px pe o linie
}

int** create_filter(char *filter_name, int *factor)
{
	int i, j, k;

	int **void_filter;
	void_filter = aloc_matrix(3, 3);

	//Crearea unui pixel care are ca valoare media tuturor punctelor din jur, inclusiv a propriei valori
	int smooth_filter[3][3] 	 = { { 1,  1,  1},
									 { 1,  1,  1},
									 { 1,  1,  1}};

	//Localizează schimbările de culoare dintr-o imagine si creeazăculori intermediare pentru a netezi marginile								 
	int blur_filter[3][3] 		 = { { 1,  2,  1},
									 { 2,  4,  2},
									 { 1,  2,  1}};

	//Detectarea diferențelor intre pixeli si accentuarea acestor diferențe							 
	int sharpen_filter[3][3] 	 = { { 0, -2,  0},
									 {-2, 11, -2},
									 { 0, -2,  0}};

	//Tot un filtru de accentuare								 
	int meanRemoval_filter[3][3] = { {-1, -1, -1},
									 {-1,  9, -1},
									 {-1, -1, -1}};

	//printf("\n\n%s %d\n\n", filter_name, strlen(filter_name));
	for(i = 0; i < 3; ++i)
		for(j = 0; j < 3; ++j)
		{
			if(strcmp(filter_name, "blur") == 0)
				void_filter[i][j] = blur_filter[i][j];

			else if(strcmp(filter_name, "smooth") == 0)
				void_filter[i][j] = smooth_filter[i][j];

			else if(strcmp(filter_name, "sharpen") == 0)
				void_filter[i][j] = sharpen_filter[i][j];

			else if(strcmp(filter_name, "mean_removal") == 0)
				void_filter[i][j] = meanRemoval_filter[i][j];
		}

	if(strcmp(filter_name, "blur") == 0)				*factor = 16;
	else if(strcmp(filter_name, "smooth") == 0)			*factor = 9;
	else if(strcmp(filter_name, "sharpen") == 0)		*factor = 3;
	else if(strcmp(filter_name, "mean_removal") == 0)	*factor = 1;
	
	return void_filter;
}

int** apply_filter(int rank, int** filter, int no_processes, int *lines_processed, int rows, int cols, int ** picture, int factor)
{
	int i, j, px;

	for(i = 1; i < rows - 1; ++i)
		for(j = 1; j < cols - 1; ++j)
		{
			px    = filter[0][0] * picture[i - 1][j - 1] +
					filter[0][1] * picture[i - 1][j    ] +
					filter[0][2] * picture[i - 1][j + 1] +
					filter[1][0] * picture[i    ][j - 1] +
					filter[1][1] * picture[i    ][j    ] +
					filter[1][2] * picture[i    ][j + 1] +
					filter[2][0] * picture[i + 1][j - 1] +
					filter[2][1] * picture[i + 1][j    ] +
					filter[2][2] * picture[i + 1][j + 1] ;
			px /= factor;
			if(px > 255) picture[i][j] = 255;
			else if(px < 0)   picture[i][j] = 0;
			else picture[i][j] = px;
		}
	
	*lines_processed += (rows-2);						 

	return picture;
}	
