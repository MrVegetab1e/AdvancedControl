#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "genetic.h"

int main (int argc, char* argv[]) {
	// Create or read a new configuration.
	Ptr_config config = (Ptr_config) malloc (sizeof(struct Chromosome_configuration));
	config->total_chrom = 1000;
	config->max_iter_best = 200;
	config->max_iter = 1000;
	config->chrom_length = 36;


	// Call genetic main method.

	genetic_main(config);
	
	return EXIT_SUCCESS;
}