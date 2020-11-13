#include <stdio.h>
#include <stdlib.h>
#include "genetic.h"


int main (int argc, char* argv[]) {
	// Create or read a new configuration.
	Ptr_config config = (Ptr_config) malloc (sizeof(struct Chromosome_configuration));
	config->total_chrom = 40;
	config->max_iter_best = 100;
	config->max_iter = 100;
	config->chrom_length = 36;
	config->mutate_percentage = 12;


	// Call genetic main method.

	genetic_main(config);
	
	return EXIT_SUCCESS;
}