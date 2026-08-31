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

void funcao_vizinhos();

#endif
