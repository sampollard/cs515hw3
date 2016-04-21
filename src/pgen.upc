#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <assert.h>

#include "packingDNAseq.h"
#include "upc_kmer_hash.h"

/** Shared Variable Declarations **/
// IDEA: Use THREADS block size (the largest possible)
// Quick reference card: http://upc.gwu.edu/downloads/quick_ref04.pdf
shared kmer_t *cur_kmer_ptr;  // Used when building contigs

int main(int argc, char *argv[]){

	/** Local Variable Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
    /* Local Variables for file IO */
    FILE *inputFile, *outputFile;
    char *input_UFX_name;
    int64_t nKmers, total_chars_to_read, cur_chars_read;
    int64_t seek_pos = 0;
    unsigned char *working_buffer;
    /* Local Variables for Graph Construction */
    hash_table_t *hashtable;   // The buckets (hashtable->table) are shared
    memory_heap_t memory_heap; // The heap is shared, the posInHeap is local.

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

    seek_pos = nKmers / THREADS;
    if (MYTHREAD % THREADS == THREADS - 1) {
        // The last thread may need to pick up the stragglers
        // This may give a lot of extra work to the last thread if THREADS is large.
        nKmers = (nKmers / THREADS) + (nKmers % THREADS);
    } else {
        nKmers /= THREADS;
    }
    printf("Thread %d reading in %d kmers\n", MYTHREAD, nKmers); // TEST
    total_chars_to_read = nKmers * LINE_SIZE;
    working_buffer = (unsigned char*) malloc(total_chars_to_read * sizeof(unsigned char));
    inputFile = fopen(input_UFX_name, "r");
    if (inputFile == NULL) {
        perror("When reading inputFile ");
    }
    // Seek into the correct offset
    fseek(inputFile, MYTHREAD*seek_pos, SEEK_SET);
    cur_chars_read = fread(working_buffer, sizeof(unsigned char), total_chars_to_read, inputFile);
    if (cur_chars_read != total_chars_to_read) {
        printf("%d read, expected %d on thread %d\n", cur_chars_read, total_chars_to_read, MYTHREAD);
        exit(1);
    }
    fclose(inputFile);

	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
    // Collectively create the hash table
    printf("Initializing hash table...\n");
    hashtable = upc_create_hash_table(nKmers, &memory_heap);
	upc_barrier;
	constrTime += gettime();
    int64_t ptr = 0;
    char left_ext, right_ext;
    start_kmer_t *startKmersList = NULL;
    printf("Creating hash table...\n");
    int startListSz = 0;
    while (ptr < cur_chars_read) {
        /* working_buffer[ptr] is the start of the current k-mer                */
        /* so current left extension is at working_buffer[ptr+KMER_LENGTH+1]    */
        /* and current right extension is at working_buffer[ptr+KMER_LENGTH+2]  */
 
        left_ext = (char) working_buffer[ptr+KMER_LENGTH+1];
        right_ext = (char) working_buffer[ptr+KMER_LENGTH+2];
 
        /* Add k-mer to hash table */
        add_kmer(hashtable, &memory_heap, &working_buffer[ptr], left_ext, right_ext);
 
        /* Create also a list with the "start" kmers: nodes with F as left (backward) extension */
        if (left_ext == 'F') {
            addKmerToStartList(&memory_heap, &startKmersList);
            startListSz++;
        }
 
        /* Move to the next k-mer in the input working_buffer */
        ptr += LINE_SIZE;
    } 
    printf("Thread %d: Start list size = %d\n", MYTHREAD, startListSz);
    
	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	upc_barrier;
    printf("Creating output file...\n");
    char filename[255];
    sprintf(filename, "pgen%d.out", MYTHREAD);
    outputFile = fopen(filename, "w");
    printf("Created output file %s\n", filename);
 
    /* Pick start nodes from the startKmersList */
    start_kmer_t *curStartNode = startKmersList; 
    kmer_t cur_kmer;
    char cur_contig[MAXIMUM_CONTIG_SIZE], unpackedKmer[KMER_LENGTH+1];
    int64_t posInContig;
    int64_t contigID = 0;
    int64_t totBases = 0;
    int startNodeList = 0;
    while (curStartNode != NULL ) {
        startNodeList++;
        /* Need to transfer the current kmer start node from shared to local */
        upc_memget(&cur_kmer, curStartNode->kmerPtr, sizeof(kmer_t));
        // Get cur_kmer to local memory
        unpackSequence((unsigned char*) &cur_kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
        /* Initialize current contig with the seed content */
        memcpy(cur_contig ,unpackedKmer, KMER_LENGTH * sizeof(char));
        posInContig = KMER_LENGTH;
        right_ext = cur_kmer.r_ext;
 
        /* Keep adding bases while not finding a terminal node */
        while (right_ext != 'F') {
           cur_contig[posInContig] = right_ext;
           posInContig++;
           /* At position cur_contig[posInContig-KMER_LENGTH] starts the last k-mer in the current contig */
            // TODO: We need to make sure the cur_kmer_ptr is a shared pointer
           cur_kmer = lookup_kmer(hashtable, (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
           right_ext = cur_kmer.r_ext;
        }
 
        /* Print the contig since we have found the corresponding terminal node */
        cur_contig[posInContig] = '\0';
        fprintf(outputFile,"%s\n", cur_contig);
        contigID++;
        totBases += strlen(cur_contig);
        /* Move to the next start node in the list */
        curStartNode = curStartNode->next;
    } 
    printf("Thread %d: From %d startNodes generated %lld contigs with %lld total bases\n",
            MYTHREAD, startNodeList, contigID, totBases);
 
    fclose(outputFile);
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
