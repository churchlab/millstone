*************************
Introduction to Millstone
*************************

What can Millstone do for me?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Millstone 0.5 currently does:

- reference-based read alignment for multiple genomes
- single nucleotide variant calling & annotation
- structural variant calling & annotation
- visualization of variants via Jbrowse
- de-novo assembly and placement of unaligned reads into contigs
- genome versioning and creation & export of new reference genomes
- variant analysis among many genomes (i.e. searching, comparison, filtering)
- design of MAGE oligos to create or revert variants

Millstone is still in active development, and there are bound to be some bugs.
Thanks for helping us find them! Please report them at our
`github repository <https://github.com/churchlab/millstone/issues>`__.

How do I get a Millstone server of my own?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently the best way to use Millstone is through Amazon AWS. Using
Amazon allows you to avoid the complexity of installing all the
dependencies from scratch on your own server, so this should be the
quickest and easiest way to deploy for most users. It requires
registering an Amazon AWS account. For projects under 50 genomes, a
suitable Amazon instance should cost less than 2 dollars per day. It is
easy to stop and start an instance when not in use. Advanced users can
also deploy their own Millstone instance locally.

We plan to write up a more complete AWS cost-guide in the future.

What do I need to get started?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You really just need two things to use Millstone:

**One or more reference genomes (Genbank or FASTA format)**: If you are
using a FASTA genome, you obviously cannot use SNPEff's variant
annotation, so we recommend Genbank if it is available. If your genome
is on Genbank, Millstone can pull the record straight from NCBI.

*Note: Millstone is meant for smaller genomes (i.e. not H. sapiens). We
use Millstone with E. coli genomes (4.6 MB) but Millstone should work
well for most microbial genomes like Saccharomyces. Try larger genomes
at your own risk.*

**Illumina HiSeq/MiSeq FASTQ Reads for one or more samples**: We've
thoroughly tested our pipeline with paired-end data, but single-end
should work as well. You need two files per sample, one for read 1,
and one for read 2. Millstone cannot (yet) split on multi-sample
barcodes or on interleaved paired-end reads, so you'll have to do that
yourself beforehand.

*Note: Extremely high-coverage samples and short fragments with
non-overlapping reads might cause difficulties. Try at your own risk.
You might consider downsampling or cleaning the reads first*