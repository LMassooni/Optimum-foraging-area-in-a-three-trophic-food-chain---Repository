/*
This file contains all the pre-defined varibles that are common to all codes and
simulations, in addition to standart C libraries and the conection between neighbors.c and the simulations codes.
*/
#ifndef COMMON
#define COMMON

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

#define FRANDOM (rand()/(RAND_MAX+1.0)) // Random variable between 0 and 1
#define L 100 // Linear side of the lattice
#define N (L*L) // Total number of sites in the lattice
#define n_especies 4 // Number of species in the simulation (empty counts as one)
#define hx 0.5 // Probability of X leaving a descendant on a site
#define gy 1.2 // Overall efficiency of predator Y
#define dy 0.15 // Death probability of Y
#define gz 1.0 // Overall efficiency of predator Y
#define dz 0.2 // Death probability of Z
#define muy 0.0 // Probability that a mutation occour on foraging area of Y 
#define muz 0.0 // Probability that a mutation occour on foraging area of Z
#define Viz_t 232 // Total number of neighbors computed for each site

extern int vizinhos[N][Viz_t]; // Array that storage the neighbors of each site (go to neighbor.c for more information)
extern int estado[N][n_especies]; // This array stores the state of the lattice. The first index corresponds to the site
// while the second is the respective species (0 = empty, 1 = X, 2 = Y, 3 = Z).
// For example, estado[12][2] = 1 indicates that in the position 12 (x = 2, y = 1) there are 1 Y individual. Only one individual is permited per site.
extern int raio_y[N]; // arrays that storage the foraging area of each individual Y in the lattice. If the site has no Y, the variable is 0.
extern int raio_z[N]; // Same as raio_y[N] but for Z predators.

extern int raios_fixos_y[34]; // List of possible foraging area values for Y species.
// As the convergence for Y can be different from Z the values of 'raios__fixos_y' and 'raios_fixos_z' could be different for better computational costs. 


extern int raios_fixos_z[34];


extern int y_length; // Length of raios_fixos_y

extern int z_length; // length of raios_fixos_z

// Vetor mapeando o limite acumulado de cada um dos 232 vizinhos
extern int n_vizinhos[Viz_t];// This array provides a way to adress new foraging areas when mutation occours. It is a list of the groups of neighbors of a given site in the vizinhos() function.
void funcao_vizinhos();
void atualizacao();

#endif
