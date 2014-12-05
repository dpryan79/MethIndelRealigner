Note: This isn't even beta quality yet. I don't recommend using it for anything serious at the moment!

MethIndelRealigner
==================
A set of programs to attempt local methylation-aware realignment around indels. This would be most useful in cases where one is performing some type of BS-seq on a human and wants to decrease incorrect methylation calls due to indel-induced misalignments. These programs work by:
 1. Finding InDels in alignments and marking their extent. A MAPQ threshold can be applied to alignments before considering whether they support the existence of an InDel.

   ![An example ROI due to an insertion](https://raw.githubusercontent.com/dpryan79/MethIndelRealigner/master/images/InDelExample.annotated.png)

   **Figure 1** Above is an example of an ROI (region of interest), where there's an insertion of CGCG versus the reference sequence. Note that some alignments are soft-clipped due to partially or totally overlapping the insertion. Since the insertion and the sequence immediately following it are the same, some reads are able to align to the reference without soft-clipping, though their 3'-most alignments are likely incorrect (they're more likely to map to the insertion).

 2. Tracking the number of alignments supporting a given InDel and using that to filter possible ROIs (regions of interest) where realignment should be performed.
 3. Looping through the input BAM/CRAM file and finding alignments that overlap an ROI.
 4. The sequence of both the reference and alignments within a 2\*Kmer window centered at the ROI (see for example [here](http://raw.githubusercontent.com/dpryan79/MethIndelRealigner/master/images/BigWindowForKmerExample.annotated.png), with K=21) are used to create a de Bruijn graph. A graph of the ROI shown above is available in [here](http://raw.githubusercontent.com/dpryan79/MethIndelRealigner/master/images/graph.pdf). The graph starts at the top and proceed through the bottom. Note that there are 6 possible paths through this graph. Blue arrows simply denote paths encountered before during depth first search.
 5. All paths (possibly without cycles, if --noCycles is used) from 2\*Kmer before the ROI to 2\*Kmer are found and alignments covering the ROI are then aligned to these paths.
 6. Each alignment from the BAM/CRAM file is assigned to the path to which it best aligns and to which the maxmimum number of other alignments best aligned (in essence, this is an expectation maximization step wherein all alignments covering an ROI are used to weight the final alignment of each other).
 7. The paths are aligned back to the reference sequence and alignments from the BAM/CRAM file are updated as needed to modify their start positions and CIGAR strings.

   ![Post-realignment](https://raw.githubusercontent.com/dpryan79/MethIndelRealigner/master/images/RealignedExample.png)

   **Figure 2** Realigned results are more probable than those prior to alignments. Note that now there are no soft-clipped bases and all of the alignments that can align to the insertion do so.

Note that while efforts are made to keep the output file sorted, there are likely edge cases where the output isn't maintained. You can alternatively pipe to `samtools sort`, though this will often degrade performance.

MethIndelRealigner is similar to BisSNP's indel realigner function, but runs in a fraction of the time (specifically, the target creator is >10x faster and the realigner is 2-3x faster). The resulting realignments tend to be similar between the two tools.

Paths are aligned to the reference sequence using a global alignment (Needleman-Wunsch) approach. BAM/CRAM alignments are realigned to paths using a global ungapped alignment approach (since there will exist a path to which an ungapped alignment can be made). This is significantly faster than regular Needleman-Wunsch.

Note that murmur3.c and murmur3.h are C implementations of MurmurHash. The C implementation is from [Peter Scott](https://github.com/PeterScott/murmur3) and MurmurHash itself is by [Austin Appleby](https://code.google.com/p/smhasher/wiki/MurmurHash3). Both of these are in the public domain.

To Do
=====
 * Ensure validity of results on a non-trivial example!!!
 * Add examples and actual documentation to the README.md
 * Using a bunch of spinlocks seems like a wasteful way to multithread. Perhaps we can chaing wake-up between functions with condition variables.
 * During graph DFS traversal, only vertices with in-degree >1 need to be tracked. This is similar to a clever memory-saving trick that minia uses. Similarly, switching to a hash would use a little more memory but end up being faster.
 * Add CRAM support.
