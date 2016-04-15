#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"

/** Shared Variable Declarations **/
// Use THREADS block size (the largest possible)
shared memory_heap_t memory_heap; // We want upc_all_alloc?

int main(int argc, char *argv[]){

	/** Local Variable Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
    /* Local Variables for reading in files */
    FILE *inputFile;
    char *input_UFX_name;
    int64_t nKmers, total_chars_to_read, cur_chars_read;
    unsigned char *working_buffer;
    // TODO: FILE *outputFile: One file or multiple?

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
    input_UFX_name = argv[1];
    /* Figure out how many lines of the file each thread is responsible for */
    nKmers = getNumKmersInUFX(input_UFX_name);

    // TEST
    if (MYTHREAD == 0) {
        printf("Total nKmers = %d\n", nKmers);
    }

    if (MYTHREAD % THREADS == THREADS - 1) {
        // The last thread may need to pick up the stragglers
        nKmers = nKmers - (THREADS-1)*(nKmers/THREADS);
    } else {
        nKmers /= THREADS;
    }
    printf("Thread %d reading in %d kmers\n", MYTHREAD, nKmers); // TEST

    total_chars_to_read = nKmers * LINE_SIZE;
    working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
    inputFile = fopen(input_UFX_name, "r");
    cur_chars_read = fread(working_buffer, sizeof(unsigned char), total_chars_to_read , inputFile);
    fclose(inputFile);

	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
    // Collectively create the hash table
    hashtable = upc_create_hash_table(nKmers, &memory_heap);
	upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
