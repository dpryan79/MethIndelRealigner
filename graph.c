#include "realigner.h"

vertex *makeVertex(char *seq, int l) {
    vertex *v = calloc(1, sizeof(vertex));
    assert(v!=NULL);
    v->seq = strdup(seq);
    v->h = hash_seq(seq, l);
    v->nchar = -1;
    return v;
}

inline void destroyVertex(vertex *v) {
    if(v->nextVertex) v->nextVertex->prevVertex = v->prevVertex;
    if(v->prevVertex) v->prevVertex->nextVertex = v->nextVertex;
    free(v->seq);
    free(v);
}

//Essentially strcmp, but for vertices
int64_t cmpVertices(vertex *a, vertex *b) {
    if(a == NULL) return 1;
    if(b == NULL) return -1;
    if(a->h == b->h) {
        if(a->nchar == b->nchar) return strcmp(a->seq, b->seq);
        return(a->nchar - b->nchar);
    }
    return ((int64_t) a->h)-((int64_t)b->h);
}

//See if a vertex is in the linked-list
//1 = yes, 0 = no
//Note that LL is modified!
int inLL(vertex **LL, vertex *v) {
    if(*LL == NULL) return 0;
    int64_t dir = cmpVertices(v, *LL);
    if(dir<0) {
        while(dir < 0) {
            if(!(*LL)->prevVertex) break;
            *LL = (*LL)->prevVertex;
            dir = cmpVertices(v, *LL);
        }
    } else if(dir > 0) {
        while(dir > 0) {
            if(!(*LL)->nextVertex) break;
            *LL = (*LL)->nextVertex;
            dir = cmpVertices(v, *LL);
        }
    }
    if(dir == 0) return 1;
    return 0;
}

int inStack(vertex **stack, vertex *v) {
    vertex *v2 = *stack; //This should be the last node
    int64_t rv;
    while(v2 != NULL) {
        rv = cmpVertices(v2, v);
        if(rv == 0) return 1;
        v2 = v2->prevVertex;
    }
    return 0;
}

//Since this was called just after inLL, LL points to the neighboring node
void addDFSVertex(vertex **LL, vertex *toAdd) {
    vertex *v = calloc(1, sizeof(vertex));
    assert(v);
    v->seq = strdup(toAdd->seq);
    assert(v->seq);
    v->h = toAdd->h;
    v->nchar = toAdd->nchar;
    int dir = cmpVertices(*LL, v);
    if(dir == 0) return; //This can only happen if we hit the same cycle more than once during DFS

    if(*LL == NULL) {
        *LL = v;
        return;
    }
    if(dir<0) {
        v->prevVertex = *LL;
        if((*LL)->nextVertex) {
            (*LL)->nextVertex->prevVertex = v;
            v->nextVertex = (*LL)->nextVertex;
        }
        (*LL)->nextVertex = v;
    } else {
        v->nextVertex = *LL;
        if((*LL)->prevVertex) {
            v->prevVertex = (*LL)->prevVertex;
            (*LL)->prevVertex->nextVertex = v;
        }
        (*LL)->prevVertex = v;
    }
    *LL=v;
}

//We'll try a recursive DFS method
void nextDFSVertex(hashTable *ht, vertex **stack, vertex *lastVertex, int k, char finalChar, vertex **visited, vertex **cycles, int32_t depth, int32_t maxDepth) {
    if(maxDepth && depth > maxDepth) return;
    vertex *v = calloc(1, sizeof(vertex));
    int i;
    char finalChars[4] = {'A', 'T', 'N', finalChar};
    assert(v);
    v->seq = strdup((*stack)->seq);
    assert(v->seq);
    v->seq = memmove((void *) v->seq, (void *) (v->seq+1), sizeof(char)*(k-1));

    for(i=0; i<4; i++) {
        v->seq[k-1] = finalChars[i];
        v->nchar = i;
        v->h = hash_seq(v->seq, k);
        if(cmpVertices(v, lastVertex) == 0) continue;
        //Can this vertex even exist?
        if(ht_sufficientCount(ht, v->seq, k)) {
            //Have we already visited this?
            if(inLL(visited, v)) {
                //Have we visited this vertex in on this path?
                if(inStack(stack, v)) {
#ifdef DEBUG
                    fprintf(stderr, "[nextDFSVertex]\t%s -> %s [color=red];\n", (*stack)->seq, v->seq); fflush(stderr);
#endif
                    inLL(cycles, *stack);
                    addDFSVertex(cycles, *stack);
#ifdef DEBUG
                } else {
                    fprintf(stderr, "[nextDFSVertex]\t%s -> %s [color=blue];\n", (*stack)->seq, v->seq); fflush(stderr);
#endif
                }
            } else { //New vertex
#ifdef DEBUG
                fprintf(stderr, "[nextDFSVertex]\t%s -> %s;\n", (*stack)->seq, v->seq); fflush(stderr);
#endif
                addDFSVertex(visited, v);
                assert(cmpVertices(v, *visited) == 0);
                //Push v onto the FIFO
                (*stack)->nextVertex = v;
                v->prevVertex = *stack;
                *stack = v;
                nextDFSVertex(ht, stack, lastVertex, k, finalChar, visited, cycles, depth+1, maxDepth);
                //Pop v back off
                v = *stack;
                *stack = (*stack)->prevVertex;
                (*stack)->nextVertex = NULL;
            }
        }
    }
    destroyVertex(v);
}

//Destroy a linked list, regardless of which member of it we're given
void destroyDFSLL(vertex *v) {
    vertex *next=v->nextVertex, *prev = NULL;
    //Destroy backwards
    while(v != NULL) {
        prev = v->prevVertex;
        destroyVertex(v);
        v = prev;
    }
    v = next;
    //Now forward
    while(v!=NULL) {
        next = v->nextVertex;
        destroyVertex(v);
        v = next;
    }
}

//Find all of the cycles in a particular subgraph with DFS
//visitedVertices: doubly-linked, ordered by hash
//cycles: doubly-linked, ordered by hash
//stack: doubly-linked FIFO
//path length
vertex * getCycles(hashTable *ht, char *startSeq, char *endSeq, int k, char finalChar, int32_t maxDepth) {
    vertex *firstVertex, *targetVertex;
    vertex *visited, *cycles = NULL;
    vertex *fifo;

    //Create the first vertex
    firstVertex = makeVertex(startSeq, k);
    firstVertex->nchar = -1;
    visited = firstVertex;
    fifo = makeVertex(startSeq, k);
    fifo->nchar = -1;
    targetVertex = makeVertex(endSeq, k);

    //Traverse
#ifdef DEBUG
    fprintf(stderr, "[nextDFSVertex]digraph %s_%s{\n", startSeq, endSeq);
#endif
    nextDFSVertex(ht, &fifo, targetVertex, k, finalChar, &visited, &cycles, k, maxDepth);
#ifdef DEBUG
    fprintf(stderr, "[nextDFSVertex]}\n");
#endif

    //Destroy the visited linked-list
    destroyDFSLL(visited);
    destroyDFSLL(fifo);
    destroyVertex(targetVertex);

    //Backup cycles to point to the first vertex
    if(cycles) while(cycles->prevVertex != NULL) cycles = cycles->prevVertex;
    return cycles;
}

//Must eventually run destroyPaths()
paths *initPaths(int m) {
    paths *p = malloc(sizeof(paths));
    p->l = 0;
    p->m = m;
    p->path = calloc(m, sizeof(char*));
    p->len = calloc(m, sizeof(int32_t));
    p->conv = 0;
    assert(p->path);
    return p;
}

void destroyPaths(paths *p) {
    int i;
    for(i=0; i<p->l; i++) {
        free(p->path[i]);
        if(p->conv && p->conv[i]) free(p->conv[i]);
    }
    free(p->path);
    if(p->conv) free(p->conv);
    free(p->len);
    free(p);
}

void appendPath(paths **paths, vertex **stack, int k, char lastChar, int32_t extraBreadth) {
    vertex *v = *stack;
    int len = 1, i;
    char *seq;

    //Have we already added too many paths?
    if(MAXBREADTH && (*paths)->l > MAXBREADTH+extraBreadth) return;

    //get the length of the path and point v to the start
    while(v->prevVertex != NULL) {
        len++;
        v = v->prevVertex;
    }
    len += k; //Otherwise, we'll be a k-mer short!
#ifdef DEBUG
    fprintf(stderr, "[appendPath] Found a path of length %i\n", len); fflush(stdout);
#endif

    //If the path length isn't at least 2k+1 in length then it isn't valid
    if(len <= 2*k+1) return;

    seq = malloc((len+1)*sizeof(char));
    assert(seq);
    seq = strcpy(seq, v->seq);
    //Add the remaining vertices
    for(i=k; i<len-1; i++) {
        v = v->nextVertex;
        seq[i] = v->seq[k-1];
    }
    seq[len-1] = lastChar;
    seq[len] = '\0';

    if((*paths)->l+1 >= (*paths)->m) {
        //Next power of 2 (32bit)
        (*paths)->m |= ((*paths)->m)>>1;
        (*paths)->m |= ((*paths)->m)>>2;
        (*paths)->m |= ((*paths)->m)>>4;
        (*paths)->m |= ((*paths)->m)>>8;
        (*paths)->m |= ((*paths)->m)>>16;
        (*paths)->m += 1;
        (*paths)->path = realloc((*paths)->path, sizeof(char*) * (*paths)->m);
        (*paths)->len = realloc((*paths)->len, sizeof(int32_t) * (*paths)->m);
        assert((*paths)->path);
        assert((*paths)->len);
    }
    (*paths)->len[(*paths)->l] = len;
    (*paths)->path[(*paths)->l++] = seq;
#ifdef DEBUG
    fprintf(stderr, "[appendPath] Appended %s\n", (*paths)->path[(*paths)->l-1]); fflush(stderr);
#endif
}

void nextBFSVertex(hashTable *ht, vertex **stack, vertex *target, int k, char finalChar, vertex **cycles, paths **paths, int32_t depth, int32_t maxDepth, int32_t minDepth, int32_t extraBreadth) {
    if(maxDepth && depth > maxDepth) return;
    vertex *v = calloc(1, sizeof(vertex));
    int i;
    char finalChars[4] = {'A', 'T', 'N', finalChar};
    assert(v);
    v->seq = strdup((*stack)->seq);
    assert(v->seq);
    v->seq = memmove((void *) v->seq, (void *) (v->seq+1), sizeof(char)*(k-1));

    for(i=0; i<4; i++) {
        if(MAXBREADTH && (*paths)->l > MAXBREADTH+extraBreadth) {
            destroyVertex(v);
            return;
        }
        v->seq[k-1] = finalChars[i];
        v->nchar = -1; //This isn't set for target
        v->h = hash_seq(v->seq, k);
        if(cmpVertices(v, target) == 0) { //Reached the target
            //Is the path long enough?
            if(depth+1 < minDepth) {
#ifdef DEBUG
                fprintf(stderr, "[nextBFSVertex] Found a path, but it was too short.\n");
#endif
                free(v->seq);
                free(v);
                return;
            }
#ifdef DEBUG
            fprintf(stderr, "[nextBFSVertex] %s -> %s (target!)\n", (*stack)->seq, v->seq); fflush(stderr);
#endif
            appendPath(paths, stack, k, finalChars[i], extraBreadth);
#ifdef DEBUG
            if((*paths)->l) fprintf(stderr, "[nextBFSVertex] Added path %s\n", (*paths)->path[(*paths)->l-1]); fflush(stderr);
#endif
            continue;
        }
        v->nchar = i;
        if(!ht_sufficientCount(ht, v->seq, k)) continue;
#ifdef DEBUG
        fprintf(stderr, "[nextBFSVertex] %s -> %s\n", (*stack)->seq, v->seq); fflush(stderr);
#endif
        if(inLL(cycles, v)) {
#ifdef DEBUG
            fprintf(stderr, "[nextBFSVertex] %s is in a cycle!\n", v->seq); fflush(stderr);
#endif
            continue;
        }

        //Push v onto the stack and recurse
        (*stack)->nextVertex = v;
        v->prevVertex = *stack;
        *stack = v;
        nextBFSVertex(ht, stack, target, k, finalChar, cycles, paths, depth+1, maxDepth, minDepth, extraBreadth);
        //Pop v back off
        v = *stack;
        *stack = (*stack)->prevVertex;
        (*stack)->nextVertex = NULL;
    }
#ifdef DEBUG
    fprintf(stderr, "[nextBFSVertex] Back-tracking\n"); fflush(stderr);
#endif
    destroyVertex(v);
}

//Return all paths between the start and end vertex in a graph, ignoring cycles
//Use a breadth-first search
paths * getPaths(hashTable *ht, char *startSeq, char *endSeq, vertex **cycles, char finalChar, int32_t maxDepth, int32_t minDepth, int32_t extraBreadth) {
    int k = strlen(startSeq);
    vertex *stack = makeVertex(startSeq, k);
    vertex *target = makeVertex(endSeq, k);
    paths *paths = initPaths(10);

#ifdef DEBUG
    fprintf(stderr, "[getPaths] Looking for paths between %"PRId32" and %"PRId32" characters long\n", minDepth, maxDepth);
#endif
    nextBFSVertex(ht, &stack, target, k, finalChar, cycles, &paths, k, maxDepth, minDepth, extraBreadth);

#ifdef DEBUG
    int i;
    fprintf(stderr, "[getPaths] Found %i paths\n", paths->l);
    for(i=0; i<paths->l; i++) fprintf(stderr, "[getPaths] path[%i] %s\n", i, paths->path[i]);
    fflush(stderr);
#endif
    destroyVertex(stack);
    destroyVertex(target);
    return paths;
}
