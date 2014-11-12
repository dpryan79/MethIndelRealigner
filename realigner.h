#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/faidx.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>

//The window size around an InDel to use for merging
#define WINDOW 17 //This should be an option, must be odd!

/*! @typedef
 @abstract This will form an element in a linked list of possible InDel Haplotypes.
 @field tid	Chromosome ID
 @field start	0-based 5' coordinate of the InDel
 @field end	0-base 3' coordinate of the InDel
 @field seq	The sequence of the Indel. This is NULL for a deletion
 @field prev	Previous node in the linked list
 @field next	next node in the linked list
*/
typedef struct InDel InDel;
struct InDel{
    int32_t tid, start, end;
    uint32_t count;
    struct InDel *prev, *next;
};

InDel *firstTargetNode, *lastTargetNode;

typedef struct{
    uint8_t *bf;
    uint64_t mask;
    uint64_t len;
} bf;

//A doubly-linked list node, or a vertex in a graph
typedef struct vertex vertex;
struct vertex {
    uint64_t h;
    char *seq;
    short nchar; //0: A, 1: T, 2: N, 3: C/G
    vertex *prevVertex, *nextVertex;
};

/*! @typedef
 @abstract A heap to hold multiple alignments (that possibly overlap an ROI)
 @field heap    The alignments
 @field l       The length of heap (i.e., number of alignments)
 @field m       The maximum number of alignments heap can hold.
 @field start   The 5' coordinate of the ROI this heap covers
 @field end     The 3' coordinate of the ROI this heap covers
*/
typedef struct {
    bam1_t **heap;
    int32_t l, m, start, end;
} alignmentHeap;

/*! @typedef
 @abstract A structure similar to kstring_t, holding all valid paths in a graph
 @field l	The number of paths
 @field m	The maximum number of paths
 @field lens	The length of each path in bases
 @field paths	The array of paths
 @field conv	The C->T/G->A sequence as int8_t[]s
 */
typedef struct {
    int32_t l, m;
    int32_t *len;
    char **path;
    int8_t **conv;
} paths;

/* Free memory allocated for a linked-list
 *
 */
void destroyTargetNodes();

/* Initialize the target nodes
 *
 * @discussion This needs to be free()d with destroyTargetNodes()
 */
void initTargetNodes();

/* Compare two nodes
 *
 * a    The first node
 * b    The second node
 *
 * returns <0 if a comes before b, >0 if b comes before a, and 0 if they overlap
 */
int TargetNodeCmp(InDel *a, InDel *b);

/* Write the target node list to a file in BED format
 * 
 * of	Output file handle
 * hdr	Header of the input BAM file
 *
 */
void writeTargets(FILE *of, bam_hdr_t *hdr);

/* Insert a node somewhere in the linked list
 *
 * node    The node to insert
 *
 * @discussion If the node matches one already in the list, then "node" is
 * freed.
 */
void insertNode(InDel *node);

/* Finds InDels in a BAM or CRAM file, adding them to the linked list
 *
 * fp    Input BAM/CRAM file
 * hdr   The header for the BAM/CRAM file
 * minMAPQ The MAPQ threshold
 *
 * discussion The linked list will need to be destroyed with destroyNodes().
 * Setting minMAPQ=-1 will use all alignments.
 */
void findInDels(htsFile *fp, bam_hdr_t *hdr, int minMAPQ);

/* Filter the linked list by node depth
 *
 * depth    The minimum read-depth of a node
 *
 * returns the number of nodes
 */
uint32_t depthFilter(int depth);

//graph.c
//Create a vertex, nchar is set to -1
//l is the length of seq
//mask is the bloom filter mask, such that subtracting two hashes doesn't overflow an int64_t
//This must be destroyed!
vertex *makeVertex(char *seq, int l, uint64_t mask);

//Given a vertex in a doubly-linked list, destroy the entire list
void destroyDFSLL(vertex *v);

//bf is the bloom filter, startSeq/endSeq are the sequences of the first/last vertices, k is the k-mer size, finalChar is C or G, if this is a G->a or C->T graph, respectively
//The output must be destroyed with destroyDFSLL()
vertex * getCycles(bf *bf, char *startSeq, char *endSeq, int k, char finalChar);

paths * getPaths(bf *bf, char *startSeq, char *endSeq, vertex **cycles, char finalChar);

void destroyPaths(paths *p);


//alignmentHeap.c
//Returns NULL on error
//size is the size of the heap
alignmentHeap * alignmentHeap_init(int size);
void alignmentHeap_destroy(alignmentHeap *heap);
void writeHeap(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap);
alignmentHeap * writeHeapUntil(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, int depth);

//bloomFilter.c
bf * bf_init(int32_t width, int kmer);
void bf_destroy(bf *bf);
void bf_reset(bf *bf);
inline void bf_add(bf *bf, uint64_t hash);
inline int bf_exists(bf *bf, uint64_t hash);
uint64_t hash_seq(char *seq, int len);

void realignHeap(alignmentHeap *heap, int k, faidx_t *fai);
