A set of programs to attempt local methylation-aware realignment around indels. This would be most useful in cases where one is performing some type of BS-seq on a human and wants to decrease incorrect methylation calls due to indel-induced misalignments.

This package uses SSW (http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0082138) to perform the actual alignments. All credit for speed should go to its authors.

Note that murmur3.c and murmur3.h are C implementations of MurmurHash. The C implementation is from {Peter Scott}[https://github.com/PeterScott/murmur3] and MurmurHash itself is by {Austin Appleby}[https://code.google.com/p/smhasher/wiki/MurmurHash3]. Both of these are in the public domain.

To Do
-----
 * Support multiple threads in the realignHeapCore.
 * Compare time requirements with BisSNP
 * Ensure validity of results on a non-trivial example!!!
 * Allow users to modify the maxmimum depth
 * Move part of realigner.c to a different file
 * Add examples and actual documentation to the README.md
 * Use a bloom filter to hold visited nodes during graph DFS traversal
