Note: This isn't even beta quality yet. I don't recommend using it for anything serious at the moment!

MethIndelRealigner
==================
A set of programs to attempt local methylation-aware realignment around indels. This would be most useful in cases where one is performing some type of BS-seq on a human and wants to decrease incorrect methylation calls due to indel-induced misalignments. These programs work by:
 1. Finding InDels in alignments and marking their extent. A MAPQ threshold can be applied to alignments before considering whether they support the existence of an InDel.
 2. Tracking the number of alignments supporting a given InDel and using that to filter possible ROIs (regions of interest) where realignment should be performed.
 3. Looping through the input BAM/CRAM file and finding alignments that overlap an ROI.
 4. The sequence of both the reference and alignments within a 2\*Kmer window centered at the ROI are used to create a de Bruijn graph.
 5. All a-cyclic paths from 2\*Kmer before the ROI to 2\*Kmer are found and alignments covering the ROI are then aligned to these paths.
 6. Each alignment from the BAM/CRAM file is assigned to the path to which it best aligns and to which the maxmimum number of other alignments best aligned (in essence, this is an expectation maximization step wherein all alignments covering an ROI are used to weight the final alignment of each other).
 7. The paths are aligned back to the reference sequence and alignments from the BAM/CRAM file are updated as needed to modify their start positions and CIGAR strings.

Paths are aligned to the reference sequence using a global alignment (Needleman-Wunsch) approach. BAM/CRAM alignments are realigned to paths using a global ungapped alignment approach (since there will exist a path to which an ungapped alignment can be made). This is significantly faster than regular Needleman-Wunsch, taking on average `O(N*(M-N)/2)` time, for alignment length N and path length M.

Note that murmur3.c and murmur3.h are C implementations of MurmurHash. The C implementation is from [Peter Scott](https://github.com/PeterScott/murmur3) and MurmurHash itself is by [Austin Appleby](https://code.google.com/p/smhasher/wiki/MurmurHash3). Both of these are in the public domain.

To Do
=====
 * Need to sort a heap before writing it to disk to ensure proper order is maintained
 * Support multiple threads in realignHeapCore:
    * findBestAlignment
    * updateAlignment
    * alignPaths2Ref
    * alignReads2Paths
 * Compare time requirements with BisSNP
 * Ensure validity of results on a non-trivial example!!!
 * Add examples and actual documentation to the README.md
 * During graph DFS traversal, only vertices with in-degree >1 need to be tracked. This is similar to a clever memory-saving trick that minia uses. Similarly, switching to a hash would use a little more memory but end up being faster.
 * Add a license (probably MIT-style)
 * Add CRAM support.
 * Determine an optimal k-mer size that generally finds the fewest paths traversing an ROI while cycles are allowed. This can then become the default.
