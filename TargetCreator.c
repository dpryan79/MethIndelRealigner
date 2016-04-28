#include "realigner.h"

//Initialize first and last
void initTargetNodes() {
    firstTargetNode = calloc(1, sizeof(InDel));
    firstTargetNode->tid = -2;
    lastTargetNode = firstTargetNode;
}

//Write the targets to a file
void writeTargets(FILE *of, bam_hdr_t *hdr) {
    lastTargetNode = firstTargetNode->next;
    while(lastTargetNode != NULL) {
        fprintf(of, "%s\t%" PRId32 "\t%" PRId32 "\n", hdr->target_name[lastTargetNode->tid], lastTargetNode->start, lastTargetNode->end);
        lastTargetNode = lastTargetNode->next;
    }
    lastTargetNode = firstTargetNode;
}

/* Compare two nodes
 *
 * a    The first node
 * b    The second node
 * k    The K-mer size
 *
 * returns <0 if a comes before b, >0 if b comes before a, and 0 if they overlap
 */
int TargetNodeCmp(InDel *a, InDel *b, int k) {
    if(a->tid == b->tid) {
        if(a->end+k < b->start-k) return -1;
        if(a->start-k > b->end+k) return 1;
        return 0;
    } else {
        return a->tid - b->tid;
    }
}

/* Free memory allocated for a node
 *
 * node   The node to free
 */
inline void destroyNode(InDel *node) {
    free(node);
}

/* Free memory allocated for a linked-list
 *
 * node    Where to start
 */
void destroyTargetNodes() {
    InDel *del;
    lastTargetNode = firstTargetNode;
    while(lastTargetNode->next != NULL) {
        del = lastTargetNode;
        lastTargetNode = lastTargetNode->next;
        destroyNode(del);
    }
    destroyNode(lastTargetNode);
}

/* Remove a node from the linked list
 *
 * node   The node to remove
 *
 * discussion Ensure that last doesn't point to this node!
 */
inline void removeNode(InDel *node) {
    node->prev->next = node->next;
    if(node->next) node->next->prev = node->prev;
    destroyNode(node);
}

/* Merge overlapping nodes
 *
 * node   The node to merge (onto last)
 * k      The K-mer size
 *
 */
void mergeNodes(InDel *node, int k) {
    InDel *nextNode = lastTargetNode->next;

    //Update last
    lastTargetNode->start = (node->start < lastTargetNode->start) ? node->start : lastTargetNode->start;
    lastTargetNode->end = (lastTargetNode->end < node->end) ? node->end : lastTargetNode->end;
    lastTargetNode->count = (lastTargetNode->count + node->count > lastTargetNode->count) ? lastTargetNode->count + node->count : lastTargetNode->count;

    if(nextNode != NULL) {
        while(TargetNodeCmp(lastTargetNode, nextNode, k) == 0) {
            lastTargetNode->end = (lastTargetNode->end < nextNode->end) ? nextNode->end : lastTargetNode->end;
            lastTargetNode->count = (lastTargetNode->count + nextNode->count > lastTargetNode->count) ? lastTargetNode->count + nextNode->count : lastTargetNode->count;
            removeNode(nextNode);
            nextNode = lastTargetNode->next;
            if(nextNode == NULL) break;
        }
    }
}

/* Insert a node into a linked list following another node
 *
 * node   The node to insert
 * target The node after which to insert
 */
void insertAfter(InDel *node, InDel *target) {
    node->prev = target;
    node->next = target->next;
    if(target->next != NULL) target->next->prev = node;
    target->next = node;
}

/* Insert a node somewhere in the linked list
 *
 * node    The node to insert
 *
 * @discussion If the node matches one already in the list, then "node" is
 * freed.
 */
void insertNode(InDel *node, int k) {
    int direction;

    //Deal with the first node
    if(lastTargetNode == firstTargetNode && firstTargetNode->next == NULL) {
        insertAfter(node, firstTargetNode);
        lastTargetNode = node;
        return;
    }

    direction = TargetNodeCmp(node, lastTargetNode, k);

    //Always come from the left
    while(direction>=0 && lastTargetNode->next != NULL) {
        lastTargetNode = lastTargetNode->next;
        direction = TargetNodeCmp(node, lastTargetNode, k);
    }
    if(direction > 0) {
        insertAfter(node, lastTargetNode);
        lastTargetNode = node;
        return;
    } else if(direction == 0) {
        mergeNodes(node, k);
        destroyNode(node);
        return;
    }

    assert(direction < 0);
    while(lastTargetNode->prev != NULL && direction < 0) {
        lastTargetNode = lastTargetNode->prev;
        direction = TargetNodeCmp(node, lastTargetNode, k);
    }
    if(direction != 0) {
        insertAfter(node, lastTargetNode);
        lastTargetNode = node;
        return;
    } else if(direction == 0) {
        mergeNodes(node, k);
        destroyNode(node);
        return;
    }
    assert(1==0);
}

/* Make a node containing an InDel
 *
 * b     The input read
 * cigar_op_num The operation number of the Insertion/Deletion
 *
 * returns a node, that must be either inserted into the linked list or
 * destroyed with destroyNode()
 */
InDel *makeNode(bam1_t *b, int cigar_op_num) {
    int i, op, oplen, quit = 0, type;
    int32_t start = b->core.pos-1;
    int32_t end;
    uint32_t *cigar = bam_get_cigar(b);
    InDel *node;

    for(i=0; i<cigar_op_num; i++) {
        oplen = bam_cigar_oplen(cigar[i]);
        type = bam_cigar_type(bam_cigar_op(cigar[i]));
        if(type & 2) start += oplen;
    }

    end = ++start;
    for(i=cigar_op_num; i<b->core.n_cigar; i++) {
        op = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        switch(op) {
        case 1: //I
        case 2: //D
            end = (end>start+oplen) ? end : start+oplen;
            break;
        default :
            quit = 1;
            break;
        }
        if(quit) break;
    }

    //Make the node
    node = calloc(1, sizeof(InDel));
    node->tid = b->core.tid;
    node->start = start;
    node->end = end;
    node->count = 1;

    return node;
}

/* Finds InDels in a BAM or CRAM file, adding them to the linked list
 *
 * fp    Input BAM/CRAM file
 * hdr   The header for the BAM/CRAM file
 * k     The K-mer size
 *
 * discussion The linked list will need to be destroyed with destroyNodes()
 */
void findInDels(htsFile *fp, bam_hdr_t *hdr, int minMAPQ, int k) {
    bam1_t *b = bam_init1();
    int i, op;
    InDel *node;
    uint32_t *cigar;

    while(sam_read1(fp, hdr, b) > 0) {
        if(b->core.qual < minMAPQ) continue;
        cigar = bam_get_cigar(b);
        for(i=0; i<b->core.n_cigar; i++) {
            op = bam_cigar_op(cigar[i]);
            if(op == 1 || op == 2) {
                node = makeNode(b, i);
                if(node == NULL) goto quit;
                insertNode(node, k);
                while(++i < b->core.n_cigar) { //Skip adjacent D/I operations
                    op = bam_cigar_op(cigar[i]);
                    if(op != 1 && op != 2) break;
                    continue;
                }
            }
        }
    }

    //Ensure that all ROIs are at least k apart
    lastTargetNode = firstTargetNode->next;
    while(lastTargetNode->next) {
        i = TargetNodeCmp(lastTargetNode,lastTargetNode->next, k);
        assert(i<=0);
        if(i==0) {
            lastTargetNode->end = lastTargetNode->next->end;
            lastTargetNode->count += (lastTargetNode->count+lastTargetNode->next->count > lastTargetNode->count)?lastTargetNode->next->count:0xFFFFFFFF;
            removeNode(lastTargetNode->next);
        } else {
            lastTargetNode = lastTargetNode->next;
        }
    }

quit:
    bam_destroy1(b);
}

/* Filter the linked list by node depth
 *
 * depth    The minimum read-depth of a node
 *
 * returns the number of nodes
 */
uint32_t depthFilter(int depth) {
    uint32_t total = 0;
    lastTargetNode = firstTargetNode->next;

    while(lastTargetNode != NULL) {
        if(lastTargetNode->count < depth) {
            lastTargetNode = lastTargetNode->prev;
            removeNode(lastTargetNode->next);
        } else {
            total++;
        }
        lastTargetNode = lastTargetNode->next;
    }
    return total;
}
