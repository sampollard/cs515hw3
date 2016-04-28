#ifndef UPC_KMER_HASH_H
#define UPC_KMER_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>

////////////////////////////////////////////
// Include contig_generation.h here because 
// we are modifying it for upc
////////////////////////////////////////////
#ifndef CONTIG_GENERATION_H
#define CONTIG_GENERATION_H

#include <upc.h>
#include <sys/stat.h>
#ifndef MAXIMUM_CONTIG_SIZE
#define MAXIMUM_CONTIG_SIZE 100000
#endif

#ifndef KMER_LENGTH
#define KMER_LENGTH 19
#endif

#ifndef LOAD_FACTOR
#define LOAD_FACTOR 1
#endif

#ifndef LINE_SIZE
#define LINE_SIZE (KMER_LENGTH+4)
#endif

static double gettime(void) {
    struct timeval tv;
    if (gettimeofday(&tv, NULL)) {
	perror("gettimeofday");
	abort();
    }
   return ((double)tv.tv_sec) + tv.tv_usec/1000000.0;
}

/* K-mer data structure */
typedef struct kmer_t kmer_t;
struct kmer_t{
   char kmer[KMER_PACKED_LENGTH];
   char l_ext;
   char r_ext;
   kmer_t *next;
};

/* Start k-mer data structure */
typedef struct start_kmer_t start_kmer_t;
struct start_kmer_t{
   shared kmer_t *kmerPtr;
   start_kmer_t *next;
};

/* Bucket data structure */
typedef struct bucket_t bucket_t;
struct bucket_t{
   shared kmer_t *head;          // Pointer to the first entry of that bucket
};

/* Hash table data structure */
typedef struct hash_table_t hash_table_t;
struct hash_table_t {
   int64_t size;            // Size of the hash table
   shared bucket_t *table;  // Entries of the hash table are pointers to buckets
};

/* Memory heap data structure */
typedef struct memory_heap_t memory_heap_t;
struct memory_heap_t {
   shared kmer_t *heap;
   int64_t posInHeap;
};

/* Returns the number of UFX kmers in a file */
int64_t getNumKmersInUFX(const char *filename) {
   FILE *f = fopen(filename, "r");
   if (f == NULL) {
      fprintf(stderr, "Could not open %s for reading!\n", filename);
      return -1;
   }
   char firstLine[ LINE_SIZE+1 ];
   firstLine[LINE_SIZE] = '\0';
   if (fread(firstLine, sizeof(char), LINE_SIZE, f) != LINE_SIZE) {
      fprintf(stderr, "Could not read %d bytes!\n", LINE_SIZE);
      return -2;
   }
   // check structure and size of kmer is correct!
   if (firstLine[LINE_SIZE] != '\0') {
      fprintf(stderr, "UFX text file is an unexpected line length for kmer length %d\n", KMER_LENGTH);
      return -3;
   }
   if (firstLine[KMER_LENGTH] != ' ' && firstLine[KMER_LENGTH] != '\t') {
      fprintf(stderr, "Unexpected format for firstLine '%s'\n", firstLine);
      return -4;
   }

   struct stat buf;
   int fd = fileno(f);
   if (fstat(fd, &buf) != 0) {
      fprintf(stderr, "Could not stat %s\n", filename);
      return -5;
   }
   int64_t totalSize = buf.st_size;
   if (totalSize % LINE_SIZE != 0) {
      fprintf(stderr, "UFX file is not a multiple of %d bytes for kmer length %d\n", LINE_SIZE, KMER_LENGTH);
      return -6;
   }
   fclose(f);
   int64_t numKmers = totalSize / LINE_SIZE;
   printf("Detected %lld kmers in text UFX file: %s\n", numKmers, filename);
   return numKmers;
}
#endif // CONTIG_GENERATION_H
////////////////////////////////////////////
// End contig_generation
////////////////////////////////////////////

/* Creates a hash table and (pre)allocates memory for the memory heap */
hash_table_t* upc_create_hash_table(int64_t nEntries, memory_heap_t *memory_heap)
{
   hash_table_t *result;
   int64_t n_buckets = nEntries * LOAD_FACTOR;

   // Each process gets its own pointer to its own chunk of the hash table
   result = (hash_table_t*) malloc(sizeof(hash_table_t));
   result->size = n_buckets;
   // There is only one table spread across all processes
   // XXX: nEntries is only enough for one thread.
   result->table = (shared bucket_t*) upc_all_alloc(n_buckets, sizeof(bucket_t));
   
   if (result->table == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the hash table: %lld buckets of %lu bytes\n", n_buckets, sizeof(bucket_t));
      exit(1);
   }
   
   // Just as with table, all processes can access the heap
   memory_heap->heap = (shared kmer_t *) upc_all_alloc(
        THREADS, ((nEntries/THREADS)+1)*sizeof(kmer_t));
   if (memory_heap->heap == NULL) {
      fprintf(stderr, "ERROR: Could not allocate memory for the heap!\n");
      exit(1);
   }
   memory_heap->posInHeap = MYTHREAD * ((nEntries / THREADS) + 1);
   
   return result;
}
/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
kmer_t lookup_kmer(hash_table_t *hashtable, const unsigned char *kmer)
{
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
   bucket_t cur_bucket;
   kmer_t result;
   static shared kmer_t *result_ptr; // hashtable->table is a shared object
   cur_bucket = hashtable->table[hashval];
   result_ptr = cur_bucket.head;
   upc_memget(&result, result_ptr, sizeof(kmer_t));
   
   /* Get the kmer from the heap and copy it to local memory */
   for (; result_ptr!=NULL; ) {
      if ( memcmp(packedKmer, result.kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
         return result;
      }
      result_ptr = result_ptr->next;
      upc_memget(&result, result_ptr, sizeof(kmer_t));
   }
   // This is never checked for in pgen.upc or serial.c
   return result;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
int add_kmer(hash_table_t *hashtable, memory_heap_t *memory_heap, const unsigned char *kmer, char left_ext, char right_ext, upc_lock_t *lock)
{
   /* Pack a k-mer sequence appropriately */
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(hashtable->size, (char*) packedKmer);
   int64_t pos = memory_heap->posInHeap;
   //printf("%d: packed = %c%c%c\n", MYTHREAD, packedKmer[0], packedKmer[1], packedKmer[2]);
   
   /* Add the contents to the appropriate kmer struct in the heap */
   /* put : private -> shared */
   upc_memput(&(memory_heap->heap[pos]).kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
   upc_memput(&((memory_heap->heap[pos]).l_ext), &left_ext, sizeof(char));
   upc_memput(&((memory_heap->heap[pos]).r_ext), &right_ext, sizeof(char));
//   printf("%d: Packed kmer = %c%c%c%c%c;%c%c;\n", MYTHREAD,
//        (memory_heap->heap[pos]).kmer[0], (memory_heap->heap[pos]).kmer[1],
//        (memory_heap->heap[pos]).kmer[2],(memory_heap->heap[pos]).kmer[3],(memory_heap->heap[pos]).kmer[4],
//        (memory_heap->heap[pos]).l_ext, (memory_heap->heap[pos]).r_ext);
   
   // TODO: Deal with the bucket stuff after the kmers have been added to a heap
   /* Fix the next pointer to point to the appropriate kmer struct */
   upc_lock(lock);
   (memory_heap->heap[pos]).next = hashtable->table[hashval].head;
   /* Fix the head pointer of the appropriate bucket to point to the current kmer */
   // The contention is at table[hashval]
   hashtable->table[hashval].head = &(memory_heap->heap[pos]);
   upc_unlock(lock);
   
   /* Increase the heap pointer */
   memory_heap->posInHeap++;
   
   return 0;
}

/* Adds a k-mer in the start list by using the memory heap (the k-mer was "just added" in the memory heap at position posInHeap - 1) */
/* The startKmersList is local, but the elements that the entries point to (kmerPtr) are shared */
void addKmerToStartList(memory_heap_t *memory_heap, start_kmer_t **startKmersList)
{
   start_kmer_t *new_entry;
   shared kmer_t *ptrToKmer;
   
   int64_t prevPosInHeap = memory_heap->posInHeap - 1;
   ptrToKmer = &(memory_heap->heap[prevPosInHeap]);
   new_entry = (start_kmer_t*) malloc(sizeof(start_kmer_t));
   new_entry->next = (*startKmersList);
   new_entry->kmerPtr = ptrToKmer; // Ptr to shared
   (*startKmersList) = new_entry;
}

/* Deallocation functions */
/* These need to be called once only */
int dealloc_heap(memory_heap_t *memory_heap)
{
   upc_free(memory_heap->heap);
   return 0;
}

int dealloc_hashtable(hash_table_t *hashtable)
{
   upc_free(hashtable->table);
   return 0;
}


#endif // KMER_HASH_H
