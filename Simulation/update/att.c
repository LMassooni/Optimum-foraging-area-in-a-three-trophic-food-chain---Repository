/*

This file contains the simulation process (update of the lattice and measurements) to an arbitrary analisys
(the dynamic is set to be the most arbitrary as possible. To reproduce some of the results you might change 
some of the lines bellow to match the respective result.)
In this file, we need to use the informations of "common.h" and neighbor.c

*/

#include "common.h"

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