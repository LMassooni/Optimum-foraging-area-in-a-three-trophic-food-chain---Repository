/*

This file contains the simulation process (update of the lattice and measurements) to an arbitrary analisys
(the dynamic is set to be the most arbitrary as possible. To reproduce some of the results you might change 
some of the lines bellow to match the respective result.)
In this file, we need to use the informations of "common.h" and neighbor.c

*/

#include "common.h"

//-------------------------------------------COMEÇO DO GNUPLOT----------------------------------------------------------------------//

int estado[N][n_especies] = {0}; // This array stores the state of the lattice. The first index corresponds to the site
// while the second is the respective species (0 = empty, 1 = X, 2 = Y, 3 = Z).
// For example, estado[12][2] = 1 indicates that in the position 12 (x = 2, y = 1) there are 1 Y individual. Only one individual is permited per site.
int raio_y[N] = {0}; // arrays that storage the foraging area of each individual Y in the lattice. If the site has no Y, the variable is 0.
int raio_z[N] = {0}; // Same as raio_y[N] but for Z predators.

int raios_fixos_y[34] = {
		0, 4, 8, 12, 20, 24, 28, 36, 44, 48, 56, 60, 68, 80, 88, 96, 100, 108, 
		112, 120, 128, 136, 144, 148, 160, 168, 176, 184, 192, 196, 212, 224, 232
		}; // List of possible foraging area values for Y species.
// As the convergence for Y can be different from Z the values of 'raios__fixos_y' and 'raios_fixos_z' could be different for better computational costs. 


int raios_fixos_z[34] = {
    0, 4, 8, 12, 20, 24, 28, 36, 44, 48, 56, 60, 68, 80, 88, 96, 100, 108, 
    112, 120, 128, 136, 144, 148, 160, 168, 176, 184, 192, 196, 212, 220, 224, 232
    };


int y_length = sizeof(raios_fixos_y)/sizeof(raios_fixos_y[0]); // Length of raios_fixos_y

int z_length = sizeof(raios_fixos_z)/sizeof(raios_fixos_z[0]); // length of raios_fixos_z

// Vetor mapeando o limite acumulado de cada um dos 232 vizinhos
int n_vizinhos[Viz_t] = { // This array provides a way to adress new foraging areas when mutation occours. It is a list of the groups of neighbors of a given site in the vizinhos() function. 
    4,4,4,4, 
    8,8,8,8, 
    12,12,12,12, 
    20,20,20,20,20,20,20,20,
    24,24,24,24,
    28,28,28,28, 
    36,36,36,36,36,36,36,36, 
    44,44,44,44,44,44,44,44, 
    48,48,48,48, 
    56,56,56,56,56,56,56,56, 
    60,60,60,60, 
    68,68,68,68,68,68,68,68,
    80,80,80,80,80,80,80,80,80,80,80,80,
    88,88,88,88,88,88,88,88,
    96,96,96,96,96,96,96,96,
    100,100,100,100, 
    108,108,108,108,108,108,108,108,
    112,112,112,112,
    120,120,120,120,120,120,120,120,
    128,128,128,128,128,128,128,128,
    136,136,136,136,136,136,136,136,
    144,144,144,144,144,144,144,144,
    148,148,148,148,
    160,160,160,160,160,160,160,160,160,160,160,160, 
    168,168,168,168,168,168,168,168,
    176,176,176,176,176,176,176,176,
    184,184,184,184,184,184,184,184,
    192,192,192,192,192,192,192,192,
    196,196,196,196,                
    212,212,212,212,212,212,212,212,212,212,212,212,212,212,212,212,
    220,220,220,220,220,220,220,220,
    224,224,224,224,                
    232,232,232,232,232,232,232,232 };


// ---------------------------------------------------- GNUPLOT ANIMATION -----------------------------------------// 
void initGnuplot(){
  printf("set autoscale fix\n");
  printf("unset xtics\n");
    printf("unset key\n");
    printf("set palette defined ( 0 \"green\", 1 \"yellow\", 2 \"red\", 3 \"blue\")\n");
  }
  
void printGnuplot(){
    FILE *arqd;
    arqd = fopen("espacial_xyz_mutacao.txt", "w+");
    int ix, iy;
    printf("plot '-' using 1:2:3 with image \n");
    for(int j=0;j<N;j++){
        ix = j%L;
        iy = j/L;
        for(int l=0;l<n_especies;l++){
            switch (estado[j][l]){
                case 0:
                    continue;
                case 1:
                    fprintf(arqd,"%d %d %d\n", ix, iy, l);
            }
        }
    }
printf("e\n");
fclose(arqd);
}

void printraio_y(){
  printf("set palette defined (0 \"white\", 4 \"green\", 8 \"yellow\", 12 \"red\", 20 \"blue\", 24 \"black\")\n");
    int ix, iy;
    printf("plot '-' using 1:2:3 with image \n");
    for(int j=0;j<N;j++){
        ix = j%L;
        iy = j/L;
		if(estado[j][2]==1)
			printf("%d %d %d\n", ix, iy, raio_y[j]);
        else
          	printf("%d %d %d\n", ix, iy, 0);
            }
printf("e\n");
}

void printraio_z(){
  printf("set palette defined (0 \"white\", 4 \"green\", 8 \"yellow\", 12 \"red\", 20 \"blue\", 24 \"black\")\n");
    int ix, iy;
    printf("plot '-' using 1:2:3 with image \n");
    for(int j=0;j<N;j++){
        ix = j%L;
        iy = j/L;
		if(estado[j][3]==1)
			printf("%d %d %d\n", ix, iy, raio_z[j]);
        else
          	printf("%d %d %d\n", ix, iy, 0);
            }
printf("e\n");
}

//-------------------------------------------END OF GNUPLOT----------------------------------------------------------------------//

// The update process of the simulation is carried out on 'atualizacao' function
void atualizacao(){
    int sitio; 
for(int h=0;h<N;h++){ // This for-loop will run 1 Monte Carlo Step (MCS) (N updates of the lattice)

	sitio = FRANDOM*N; // sitio is the site that will be selected randomly to be update 
	
    for(int i = 0; i<1;i++){// This for-loop exists solely so that the "break;" statements exit it, rather than the outer loop. This ensures that
    // whenever the code reaches a break, it terminates the dynamics for that site, exits the i-loop, and returns to the h-loop to select another site and re-run the dynamics.
    
// Step 1 - Populate empty sites
        if(estado[sitio][0] == 1){
        
            for(int m = 0;m<4;m++){ // Only the four nearest neighbors can leave a descendant on the central site
                if(estado[vizinhos[sitio][m]][1]!=0 && FRANDOM < hx){
                // For X to appear at the site, a neighbor must be X. 'estado[vizinhos[sitio][m]] != 0' takes care of identifying if the neighbors of "sitio" has an X specie (estado[neighbor][1] = 1 ). 
                //If this occurs and the random number is < hx, the code enters this if-block and reproduces.
                    estado[sitio][0]=0;
                    estado[sitio][1]=1;
                    break;
                }
            }
			// If there are no Xs among the neighbors of 'sitio', or if the random numbers turn out not to be < hx, the code will never enter the 'if' block. If this happens, we need
            // to terminate the dynamics for this site—hence the 'break' below.
            break;
        }
        
// Step 2: Probability of Y preying on X

		if(estado[sitio][1]==1){

		bool comido = false; // Boolean variable that will control if the predation was successful
		
            for(int b = 1; b<y_length;b++){ //Loop that will run through all the possible neighbors (raios_fixos_y)
            	int inicial = raios_fixos_y[b-1];
            	int final = raios_fixos_y[b];	 
            	int tam = final-inicial; // Number of neighbors in the same distance of the central site
            	int offset = FRANDOM*tam; // Random neighbor in the distance

				for(int s = 0;s<tam;s++){
		            	int k = inicial + ((offset + s) % tam);// Accesses neighbors at the same distance in a random circular order
		            	            // Again, for the prey at 'sitio' to be preyed upon, a neighbor must be of type Y. After checking this, we check the predation probability
		            	            // and also whether that neighboring predator has a predation radius sufficient to reach the prey. If, for example, k = 12
		            	            // but the neighbor's predation radius is 4, the final expression will evaluate to false.
					if(estado[vizinhos[sitio][k]][2]==1){
	                if(FRANDOM < gy/(double)raio_y[vizinhos[sitio][k]] && raio_y[vizinhos[sitio][k]]>=n_vizinhos[k]){
		// If all conditions are met, we enter this if-block and perform the reproduction process using the mutation radius probability. Otherwise, we move on to the next neighbor.
			            estado[sitio][1]=0;
	                    estado[sitio][2]=1;
	                    raio_y[sitio] = raio_y[vizinhos[sitio][k]]; // Atribution of the same foraging area of the ancestor

	                    comido = true; // Bolean variable changes to True
	                    
	                    double raio = FRANDOM; // This random number will be used to determine whether to increase or decrease the radius.
						if(FRANDOM<muy){ // If execution reaches here, a mutation has occurred

						// The mutation of the 1st neighbors and the 'last' neighbors needs to be handled separately
							if(raio_y[sitio]==n_vizinhos[0]){  // if the area is 4
								raio_y[sitio] = n_vizinhos[0]+4;
								break;
							}
							if(raio_y[sitio]==raios_fixos_y[y_length-1]){ // if the area is the higher possible 
								raio_y[sitio] -= 4;// This value must be change manually. For example, if the last value is 232, it must be -8, if is 224, it must be -4
								break;							
							}
							if(raio_y[sitio]== 80){ // 80 and 160 neighbors is special because it differs in 12 neighbors from the previous group (68, 148)
								if(raio>0.5){
									raio_y[sitio] = raio_y[sitio]+8; 
									break;
								}
								else{
									raio_y[sitio] = raio_y[sitio]-12;
									break;
								}
							}
							if(raio_y[sitio]== 212){ // 212  neighbors is special because it differs in 16 neighbors from the previous group (196)
								if(raio>0.5){
									raio_y[sitio] = raio_y[sitio]+8; 
									break;
								}
								else{
									raio_y[sitio] = raio_y[sitio]-16;
									break;
								}
							}
							// The rest of the groups has no special features and enter here
							if(raio>0.5){
									raio_y[sitio] = n_vizinhos[raio_y[sitio]];
									break;
								}
							else{
								if(n_vizinhos[raio_y[sitio]-1]==n_vizinhos[raio_y[sitio]-5])
									raio_y[sitio]=n_vizinhos[raio_y[sitio]-9];
								else
									raio_y[sitio]=n_vizinhos[raio_y[sitio]-5];
								break;
								}				
							}
							break;
							// All the breaks above will exit only the "same group" loop, but if the prey was sucessfull eaten, we cannot permit that the outer loop (variable 's') continue.
							// To account this, we check the bolean variable "comido". It only changes to "True" if it entries the all the conditionals of predation.
	                	}
	                	} 	
	                }
	                // If the preying was not sucessfull, comido is still false, and this loop will continue until all the neighbors were swept away or if comido == True
	                if(comido == true)
	                	break;
	            }
	      break;
        }
            
 //Step 3: death of predator Y and predation of Y by Z
        if(estado[sitio][2]==1){
        // First, we check if Y dies naturally
            if(FRANDOM < dy){ 
                estado[sitio][2]=0;
                estado[sitio][0]=1;
                raio_y[sitio] = 0;
                break;
            }
            
            else{
            // If it doesn't die, we enter this else block. The dynamics here are exactly the same as the predation of X by Y.

            bool comido = false;
            for(int b = 1; b<z_length;b++){
            	int inicial = raios_fixos_z[b-1];
            	int final = raios_fixos_z[b]; 
            	int tam = final-inicial;
            	int offset = rand() % tam;
				for(int s = 0;s<tam;s++){
	            	int k = inicial + ((offset + s) % tam);
	            	if(estado[vizinhos[sitio][k]][3]==1){
		                if(FRANDOM < gz/(double)raio_z[vizinhos[sitio][k]] && raio_z[vizinhos[sitio][k]]>=n_vizinhos[k]){
		                	estado[sitio][2]=0;
		                    estado[sitio][3]=1;
		                    raio_y[sitio] = 0;
		                    comido = true;
		                    raio_z[sitio] = raio_z[vizinhos[sitio][k]];
		                    double raio = FRANDOM;
							if(FRANDOM<muz){
								if(raio_z[sitio]==n_vizinhos[0]){
									raio_z[sitio] = n_vizinhos[0]+4;
									break;
									}
								if(raio_z[sitio]==n_vizinhos[Viz_t-1]){
									raio_z[sitio] -= 8;
									break;							
								}
								if(raio_z[sitio]==80 || raio_z[sitio] == 160){
									if(raio>0.5){
										raio_z[sitio] = raio_z[sitio]+8;
										break; 
									}
									else{
										raio_z[sitio] = raio_z[sitio]-12;
										break;
									}
								}
								if(raio_z[sitio]==212){
									if(raio>0.5){
										raio_z[sitio] = raio_z[sitio]+8;
										break; 
									}
									else{
										raio_z[sitio] = raio_z[sitio]-16;
										break;
									}
								}

								if(raio>0.5){
		        					raio_z[sitio] = n_vizinhos[raio_z[sitio]];
		        					break;
		        				}
						        else{
						        	if(n_vizinhos[raio_z[sitio]-1]==n_vizinhos[raio_z[sitio]-5])
						        		raio_z[sitio]=n_vizinhos[raio_z[sitio]-9];
						        	else
						        		raio_z[sitio]=n_vizinhos[raio_z[sitio]-5];
						        	break;
						        	}				
						        }
						        break;
		                	}
	                	} 	
	                }
	          if(comido == true)
	          	break;  
             }
           break;
         }
     }

// Finally, if the site "sitio" has a Z species, it can only die.
        if(estado[sitio][3]==1 && FRANDOM<dz){ 
            estado[sitio][3]=0;
            estado[sitio][0]=1;
            raio_z[sitio] = 0;
            break;
			}
        }
    }
}



// This function calculates simple measures, such as the average populations and foraging areas

void medidas(double media[3], double raios[2]){

	media[0] = media[1] = media[2] = raios[0] = raios[1] = 0;

	for(int i=0;i<N;i++){

// Check if the unity carrying capacity is working
	if((estado[i][1] > 1 || estado[i][2] > 1 || estado[i][3] > 1) ||
	 (estado[i][1] == 1 && estado[i][2] == 1) ||
	 (estado[i][1] == 1 && estado[i][3] == 1) || 
	 (estado[i][2] == 1 && estado[i][3] == 1) ){
		printf("ERROR! - More than 1 individual per site (site - %d)\n", i);
	}
		media[0] += estado[i][1];

		media[1] += estado[i][2];

		media[2] += estado[i][3];

		raios[0] += raio_y[i];

		raios[1] += raio_z[i];

	}
	raios[0] = raios[0]/media[1];
	raios[1] = raios[1]/media[2];
	media[0] = media[0]/N;
	media[1] = media[1]/N;
	media[2] = media[2]/N;
}

//This function calculates more complicated measurements - such as the frequency of each foraging area in a given population
void dispersion(double media[3], double raios[2], double dispersao_y[34], double dispersao_z[34]){

media[0] = media[1] = media[2] = raios[0] = raios[1] = 0;
int valor_raio;
for(int b = 1; b<34;b++){
	dispersao_y[b] = 0;	
	dispersao_z[b] = 0;
}
	for(int i=0;i<N;i++){
	valor_raio = 0;
	if((estado[i][1] > 1 || estado[i][2] > 1 || estado[i][3] > 1) ||
	 (estado[i][1] == 1 && estado[i][2] == 1) ||
	 (estado[i][1] == 1 && estado[i][3] == 1) || 
	 (estado[i][2] == 1 && estado[i][3] == 1) ){
		printf("ERROR! - More than 1 individual per site (site - %d)\n", i);
	}
		media[0] += estado[i][1];
		media[1] += estado[i][2];
		media[2] += estado[i][3];

//This block addresses the quantity of Y individuals with each area
		if(estado[i][2]==1){
			for(int j = 1; j<y_length;j++){
				if(raio_y[i]==raios_fixos_z[j]){
					valor_raio = j;
					break;
				}
			}
			dispersao_y[valor_raio] += 1; // array that stores the population per area for Y
		}
		
//This block addresses the quantity of Z individuals with each area
	if(estado[i][3]==1){
		for(int j = 1; j<z_length;j++){
			if(raio_z[i]==raios_fixos_z[j]){
				valor_raio = j;
				break;
			}
		}
		dispersao_z[valor_raio] += 1; // array that stores the population per area for Z

	}		
	raios[0] += raio_y[i];
	raios[1] += raio_z[i];	
	}
	for(int b = 1; b<34;b++){
		dispersao_z[b] = dispersao_z[b]/media[2];
		dispersao_y[b] = dispersao_y[b]/media[1];
		}
	raios[0] = raios[0]/media[1];
	raios[1] = raios[1]/media[2];
	media[0] = media[0]/N;
	media[1] = media[1]/N;
	media[2] = media[2]/N;
}



// Main function that will run the update and perform the desired measurements.
int main(){
	funcao_vizinhos(); // Calls the function that address the neighbors in "neighbors.c" 
    srand(time(NULL)); // Defines the seed of the random numbers according to the time of the computer

// Criation of an external file that will storage the result of the simulatiotn
    FILE *arq;
    arq = fopen("name-of-result-file.txt", "w+");

    initGnuplot(); // Function that initiates the gnuplot for animation
	double p = 10000; // Number of MCS for the measurement process
	int m = 5000; // Number of MCS of "thermalization"
	int initial_y_area = 4;
	int initial_z_area = 4;
	double media[3] = {0}; // Array that will storage the populations density (0 = X, 1 = Y, 2 = Z)
	double raios[2] = {0}; // Mean foraging area of the populations (0 = Y, 1 = Z)

 
// Initial condition
    for(int i=0; i<N;i++){
    	estado[i][0] = estado[i][1] = estado[i][2] = estado[i][3] = 0;
    	raio_z[i] = raio_y[i] = 0;
        int m = FRANDOM*4; // Selects randomly a number between 0 and 3 (for the species)
        estado[i][m] = 1;
        if(m==2){
        	raio_y[i]=initial_y_area;
        }
        if(m==3){
        	raio_z[i]=initial_z_area;
        }
    }

//Transient period. Involves in
/*
*/
// If you want to skip this step (for results such as Figure 4 and 5 that includes the time series from the initial condition)
	for(int t = 0; t<m;t++){
		atualizacao();
	}


//Simplified Update process
    for(int t=0; t<(int)p; t++){
        atualizacao();
	    medidas(media,raios);
	    fprintf(arq, "%d %lf %lf %lf %lf %lf\n",t ,media[0],media[1],media[2],raios[0],raios[1]); // Copies the variables values to the file created earlier
//If you want to see the dynamics, you can use the functions printraio_z, printraio_y and printGnuplot to see the lattice changing in real time.
//printraio_z shows the Z foraging area distribution
//printraio_y shows the Y foraging area distribution
//printGnuplot shows the species distribution.
//Only one of them will work per simulation, so remember to comment with "//" the undesired outputs

//Remove "/*" "*/" to see the real time simulation. More information on how to use gnuplot in this context will be in the ReadMe file 
/*
	    if(t%1000 == 0){
	    	printraio_z();
	    	printraio_y();
	    	printGnuplot();
	    }
*/
  	}
  	
fclose(arq);

// For more complicated dynamics, with repetition and different measures, you must use the next update process	
    FILE *file1;
    file1 = fopen("file_1.txt", "a");
    FILE *file2;
    file2 = fopen("file_2.txt", "a");
    FILE *file3;
    file3 = fopen("file_3.txt", "a");

//	double media[3] = {0};
	//Arrays that will storage the frequencies of each foraging area for Y and Z
	double dispersao_y[34] = {0};
	double dispersao_z[34] = {0};
//	double raios[2] = {0};

// Auxiliary variables
	double areay=0;
	double areaz=0;
	double areaf1=0;
	double areaf2=0;
	double dispersao_t1[34] = {0};
	double dispersao_f1[34] = {0};
	double dispersao_t2[34] = {0};
	double dispersao_f2[34] = {0};

//Number of histories (equal simulations)
	double historias = 10.0;

	for(int u = 0; u < (int)historias; u++){
		for(int i=0; i<N;i++){

			estado[i][0]=estado[i][1]=estado[i][2]=estado[i][3]=0;
			raio_y[i]=raio_z[i]=0;
		    int m = FRANDOM*4;
		    estado[i][m] = 1;
		    if(m==2){
		    	raio_y[i]=initial_y_area;
		    }
		    if(m==3){
		    	raio_z[i]=initial_z_area;
		    }
		}
		for(int t = 0; t<m;t++){
			atualizacao();
		}
		
	    for(int t=0; t<(int)p; t++){
	        atualizacao();
	        dispersion(media,raios, dispersao_y, dispersao_z);
	        areay += raios[0]; // Sums the mean Y foraging area
	        areaz += raios[1];// Sums the mean Z foraging area
	        
	        for(int b = 1;b<y_length;b++){
	        	dispersao_t1[b] += dispersao_y[b]; // Sums the density of each foraging area for Y
	        }
	        for(int b = 1;b<z_length;b++){
	        	dispersao_t2[b] += dispersao_z[b]; // Sums the density of each foraging area for Z
	        }
	  	}
	  	// Takes the average of 'p' MC
	  	areaf1 += areay/p;
	  	areaf2 += areaz/p;
	  	areay = 0;
	  	areaz = 0;
	  	for(int b = 1;b<y_length;b++){
	  		dispersao_f1[b] += dispersao_t1[b]/p;
	  		dispersao_t1[b] = 0;
		}
		for(int b = 1; b< z_length; b++){
	  		dispersao_f2[b] += dispersao_t2[b]/p;
	  		dispersao_t2[b] = 0;
	  	}
	  	//end of the history
	}
// Saves the average (of the histories) of the desired measures  	  	
	  	fprintf(file1, "%lf %lf\n",areaf1/historias ,areaf2/historias);
	  	areaf1 = areaf2 = 0;
			for(int b = 1;b<y_length;b++){
				fprintf(file2, "%d %lf\n", raios_fixos_y[b], dispersao_f1[b]/historias);
				dispersao_f1[b] = 0;
			}
			for(int b = 1;b<z_length;b++){
				fprintf(file3, "%d %lf\n", raios_fixos_z[b], dispersao_f2[b]/historias);
				dispersao_f2[b] = 0;
			}
// The block bellows works for the case when loops are perfomed beyond the histories loops.
// For example, Figure 3, 4-c and 7 used a loop in fixed foraging areas. In this case, we can perform a number of histories for the same foraging area, but
// varies the foraging area itself. In this case, it is preferable to save the results at each foraging area value, in case something goes wrong with the computer
// performing the simulation or another external factor, since this simulations might take a while to finish
/*
fclose(file1);
fclose(file2);
fclose(file3);
file1 = fopen("file1.txt", "a");
file2 = fopen("file2.txt", "a");
file3 = fopen("file3.txt", "a");
*/
fclose(file1);
fclose(file2);
fclose(file3);

//To perform some of the results shown in the paper, go to the ReadMe file for more information

return 0;
}
