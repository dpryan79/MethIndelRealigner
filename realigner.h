#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/faidx.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#define REALIGN_H 1

/*! @typedef structure of the alignment result
 @field score1      the best alignment score
 @field ref_begin1  0-based best alignment beginning position on reference
 @field ref_end1    0-based best alignment ending position on reference
 @field read_begin1 0-based best alignment beginning position on read
 @field read_end1   0-based best alignment ending position on read
 @field cigar       best alignment cigar, as in htslib
 @field cigarLen    length of the cigar string; cigarLen = 0 when the best
                    alignment path is not available
 @discussion This is basically a subset of the s_align struct in SSW, which used
 to be used for Smith-Waterman realignment to paths.
*/
typedef struct {
    uint16_t score1;
    int32_t ref_begin1;
    int32_t ref_end1;
    int32_t read_begin1;
    int32_t read_end1;
    uint32_t* cigar;
    int32_t cigarLen;
} s_align;

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

typedef struct hashTableEntry hashTableEntry;
struct hashTableEntry {
    char *seq;
    uint8_t cnt;
    struct hashTableEntry *next;
};
    
typedef struct {
    int threshold;
    uint64_t n;
    hashTableEntry **entries;
} hashTable;

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

int32_t MAXBREADTH; //The maximum number of C->T or G->A paths allowed. If more than this are found, no realignment occurs.
int32_t MAXINSERT; //The maximum insert size that will be looked for. This limits the depth of graph DFS
int32_t MINMAPQ; //Minimum MAPQ score an alignment must have to process it.
faidx_t *GLOBAL_FAI; //In order to ensure that a realignment is best, a deeply buried function needs access to the full reference sequence
int32_t NOCYCLES; //Whether to allow cycles in a path or not (default: allowed).

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
 * k    The K-mer size
 *
 * returns <0 if a comes before b, >0 if b comes before a, and 0 if they overlap
 */
int TargetNodeCmp(InDel *a, InDel *b, int k);

/* Write the target node list to a file in BED format
 * 
 * of	Output file handle
 * hdr	Header of the input BAM file
 *
 */
void writeTargets(FILE *of, bam_hdr_t *hdr);

/* Insert a node into a linked list following another node
 *
 * node   The node to insert
 * target The node after which to insert
 */
void insertAfter(InDel *node, InDel *target);

/* Insert a node somewhere in the linked list
 *
 * node The node to insert
 * k    The K-mer size
 *
 * @discussion If the node matches one already in the list, then "node" is
 * freed.
 */
void insertNode(InDel *node, int k);

/* Finds InDels in a BAM or CRAM file, adding them to the linked list
 *
 * fp    Input BAM/CRAM file
 * hdr   The header for the BAM/CRAM file
 * minMAPQ The MAPQ threshold
 * k     The k-mer size
 *
 * discussion The linked list will need to be destroyed with destroyNodes().
 * Setting minMAPQ=-1 will use all alignments.
 */
void findInDels(htsFile *fp, bam_hdr_t *hdr, int minMAPQ, int k);

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
//This must be destroyed!
vertex *makeVertex(char *seq, int l);

//Given a vertex in a doubly-linked list, destroy the entire list
void destroyDFSLL(vertex *v);

//ht is the hash table, startSeq/endSeq are the sequences of the first/last vertices, k is the k-mer size, finalChar is C or G, if this is a G->a or C->T graph, respectively
//The output must be destroyed with destroyDFSLL()
vertex * getCycles(hashTable *ht, char *startSeq, char *endSeq, int k, char finalChar, int32_t maxDepth);

paths * getPaths(hashTable *ht, char *startSeq, char *endSeq, vertex **cycles, char finalChar, int32_t maxDepth, int32_t minDepth, int32_t extraBreadth);

void destroyPaths(paths *p);


//alignmentHeap.c
//Returns NULL on error
//size is the size of the heap
alignmentHeap * alignmentHeap_init(int size);
void alignmentHeap_destroy(alignmentHeap *heap);
void writeHeap(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, int doSort);
alignmentHeap * writeHeapUntil(samFile *of, bam_hdr_t *hdr, alignmentHeap *heap, bam1_t *curb, samFile *fp);

//hashTable.c
hashTable *ht_init(int32_t width, int kmer, int threshold);
void ht_addIncrement(hashTable *ht, char *seq, int len);
void ht_addIncrementMax(hashTable *ht, char *seq, int len);
int ht_sufficientCount(hashTable *ht, char *seq, int len);
uint64_t ht_numEntries(hashTable *ht);
void ht_destroy(hashTable *ht);
uint64_t hash_seq(char *seq, int len);

//realign.c
void realignHeap(alignmentHeap *heap, int k, faidx_t *fai, int nt, int threshold);

//needlemanWunsch.c
s_align * GlobalAlignment(int8_t *ref, int32_t refLen, int8_t *path, int32_t pathLen, int k, int32_t likelyStartpos);
s_align * SemiGlobalAlignment(int8_t *ref, int32_t refLen, int8_t *path, int32_t pathLen, int32_t k);
