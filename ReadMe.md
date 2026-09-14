# Welcome to the repository of the work - Optimum foraging area in a three-trophic food chain


## This file contains the overall information to perform the simulations and analytical dynamics of the model.

# Simulation

## The simulation depends, in general, of four files:
### - common.hpp has all the common constants of the model, such as the probabilities of predation/leave a descendant, or the death probabilities. It also connects the other two files and loads the essential packages of C language. To use this file, you will need to
```
#include "common.h"
```
### in the beginning of the files.

### - neighbors.c has the process of declaring and assign the correct neighbors of each site. The highest number of neighbors used was 232. This means that each site of the simulation has 232 neighbors assigned in the array "vizinhos". If you want to add more neighbors, make sure to add the coordinates distance in the array "lis", with the x and y distance in a tuple form. The distance needs to be ordered with the same distance from the center. If you add new neighbors, remember to change the value "Viz_t" on the "common.h" file to the new total value of neighbors. This two variables (Viz_t and lis) are the only thing you will need to change.

### - att.c has the update process of the model. Details on how the code works are commented on the file 

### - main.c is the file that call the update function from att.c and the neighbors array from neighbors.c and make the desired measurements.

### (// and /* */ blocks signals comments).

## How to reproduce the simulation results of the work

### att.c is in the most general form, but to reproduce the results of the work, some changes may be needeed.
### For example, Figure 2 uses a loop on the g_z variable. To perform it, one must comment or exclude the line
```
#define gz 1.0
```  
### from the common.h file and add a loop on line 450 or 523 that runs different values of g_z. For example, one could make:
```
double gz = 0.0 // At the beginning of the code
//line 450
for(int v = 0 ; v<50 ; v++){
	//block of code that is already there
//line 494
gz = gz + 0.08
}
```
### With the only change being in the save file, where you save only the desired variable (the population density in this case) and in this case eh average over some MCS.

### Figure 3 and 4-c is made similarly, with looping over the foraging area of both Y and Z.

## REMEMBER OF ALWAYS CHECK THE VALUES OF MUTATING RATE (muz and muy)

### For Figure 5, one must save the averages foraging areas as function of time, as in the second part of att.c. The subfigures (b) and (c ) are made saving the frequencies as functions of time as well.

### Finally, Figure 7 is performed with the following block of code:

```
	int invadiu = 0; // variable that will count the number of sucessfull invasions for a given resident and invader population
	double pop_inicial, pop_invasor, pop_x, pop_z, pop_inicial_h, pop_invasor_h, pop_x_h, pop_z_h=0;
	int raio_inicial;	
for(int x = 4; x<7;x++){
	raio_inicial = raios_fixos[x]; // foraging area of the resisdent population
	for(int u = 1; u<tamanho; u++){
		int raio_invasor = raios_fixos[u]; // foraging area of invader population
//....
// (PART OF INITIAL CONDITION AND LOOP IN HISTORIES)
//....
			for(int k = 0;k<m;k++){ // m = 5000 in this case
				atualizacao();
			}
			medidas(media, raios, raio_inicial, raio_invasor);

			//This block takes 1% of the Y population and changes the foraging area to the invader value
			for(int j = 0; j<(int)(media[1]*N/100);j++){
				for(int h = 0; h<N;h++){
					int sitio = FRANDOM*N;
					if(estado[sitio][2]==1){
						raio_y[sitio] = raio_invasor;
						break;
					}
				}

			}

			// Returns the update
			medidas(media, raios, raio_inicial, raio_invasor);
			for(int b = 0;b<(int)m;b++){
				atualizacao();
				medidas(media, raios, raio_inicial, raio_invasor);
				pop_inicial += raios[0];
				pop_invasor += raios[1];
				pop_x += media[0];
				pop_z += media[2];
			}
			pop_inicial_h += pop_inicial/m;
			pop_invasor_h += pop_invasor/m;
			pop_x_h += pop_x/m;
			pop_z_h += pop_z/m;
			// The condition to a sucessfull invasion is that the  invader population is greater than the resident population, and is
			// greater than a treshold of 0.05
			if(pop_invasor > pop_inicial && pop_invasor > 0.05){
				invadiu +=1;
			}
			pop_inicial = 0;
			pop_invasor = 0;
			pop_x = 0;
			pop_z = 0;

		 	}
		fprintf(arq1, "%d %d %d\n", raios_fixos[x],raios_fixos[u], invadiu);
		invadiu = 0;
		pop_inicial_h = 0;
		pop_invasor_h = 0;
		pop_x_h = 0;
		pop_z_h = 0;
	}
``` 
## To run any C file, in the linux terminal you must execute:
```
gcc main.c update/*.c -lm -O3 -o %name of executable%

./name of executable
```


# Analytical approach
