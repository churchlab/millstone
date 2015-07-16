********
Variants
********

All sample, reference genome, and alignments are listed in the *Data*
view (the toggle switch in the top left). Clicking over to the *Analyze*
view will allow you to filter through multi-sample variants and view
their aligned reads. Use the dropdowns on the left below the
Data/Analyze toggle to select your alignment and your reference genome
and choose *Variants*. Once the alignments are complete, you should see
a list of all variants that have been identified across all samples.

Cast vs. Melted
---------------

There are two ways to view variants.

**Cast**: Cast displays a summary row for one variant across all
samples. You can see how many samples the variant is present in, as well
as the variant's effects.

**Melted**: 'Melting' the view shows one row for every combination of
sample and and variant. It essentially multiplies the rows by the number
of samples, so you can see data specific to individual samples. If a
variant is not called in a sample, it's *Alt* column will be blank.

Links
-----

There are three link icons next to every sample.

-  *The magnifying glass* icon 'zooms in' to the melted view for that
   variant across all samples.
-  *The read alignment* icon shows how individual fastq reads align
   around a variant. It is useful for doing visual QC on an alignment,
   to make sure your reads are properly aligned around your variant.
-  *The bar graph* icon shows the coverage of your reads. Areas of high
   or low coverage might be of interest, and this view is more compact,
   which makes it easier to compare multiple samples.

*Note: If an icon is gray in the Cast view* it is disabled because it is
too intensive to display many samples simultaneously. Zoom into the
variant (with the magnifying glass) and inspect individual samples. You
can manually add and remove tracks in Jbrowse via the track list on the
left.

More information about using JBrowse and understanding its visualization
can be found at its `website <http://jbrowse.org/>`__.

Fields and Filtering
--------------------

Millstone uses a simple language to understand query syntax for
filtering variants.

*Note: Currently some of the field names can be confusing.* A list of
all available fields can be found with the *Fields...* button. The
default column names don't always correspond to the internal field
names. There isn't currently a well-documented list of what each field
means, but most of them are documented in the `VCF
specification <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`__.
The ``INFO_EFF_*`` fields come from
`SnpEFF <http://snpeff.sourceforge.net/SnpEff_manual.html#output>`__.

Examples
^^^^^^^^

If you want to look at all variants in a certain gene:

``INFO_EFF_GENE = tolC``

If you want to look at all variants that have strong or moderate
predicted phenotypic effects:

``INFO_EFF_IMPACT = HIGH | INFO_EFF_IMPACT = LOW``

If you want to look in a certain region:

``CHROM = NC_000913 & POSITION > 500 & POSITION < 1000``

Marginal Calls
--------------

We always run variant calling as diploid, even for haploid organisms
like E. coli, so that some poorly-supported variants appear
heterozygous. This allows marginal calls to be made in cases where only
a portion of the reads show a SNV, in cases of regional duplications or
if reads map to a non-unique region of the genome. Such marginal calls
have an orange fraction icon in their ALT column, and can also be
filtered on by using:

``IS_HET = TRUE`` or ``IS_HET = FALSE``

Additionally, the ``GT_TYPE`` field is another way to distinguish
marginal from strongly called variants. ``GT_TYPE`` can take values
between 0 or 1 for each sample/variant combination:

-  0 means the variant was called as reference in the sample
-  1 means the variant was called as heterozygous (i.e. marginal) in the
   sample
-  2 means the variant was called as homozygous (well-supported) in the
   sample

If you'd like to filter on only well-supported variants that have
moderate to strong affects on genes, you can use the filter:

``GT_TYPE = 2 & (INFO_EFF_IMPACT = HIGH | INFO_EFF_IMPACT = MODERATE)``

Variant Sets
------------

Variant sets are a way to group variants after filtering. The sets
created by default correspond to regions where the alignment had
problems; either there was insufficient coverage, no coverage, too much
coverage, or poor mapping quality (corresponding perhaps to regions that
are non-unique).

You can also create your own sets to group interesting variants, or
those whose alignments you'd like to examine by eye.

Creating a blank set
^^^^^^^^^^^^^^^^^^^^

You can create your own blank sets from the Sets tab in the *Data* view.
After creating a set, you can add variants to it in the *Analyze* view
using the checkboxes and the master checkbox dropdown on the left.

Uploading a set from a VCF file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also upload a variant set from a VCF file. Only the first 5
columns of the VCF will be used. The file must be tab delimited. Here is
an example:

::

    #CHROM          POS ID  REF ALT
    NC_000913   2242    .   G   A
    NC_000913   76  .   C   A
    NC_000913   3170    .   T   C
    NC_000913   1623    .   G   C
    NC_000913   3879    .   A   G
    NC_000913   3112    .   A   T
    NC_000913   1577    .   C   T
    NC_000913   5352    .   G   A
    NC_000913   4386    .   A   T
    NC_000913   1167    .   G   T
    NC_000913   5425    .   T   A
    NC_000913   951 .   C   A
    NC_000913   3993    .   A   G
    NC_000913   226 .   G   C
    NC_000913   2939    .   T   G
    NC_000913   92  .   C   A
    NC_000913   5563    .   A   C
    NC_000913   4446    .   A   C
    NC_000913   607 .   A   G
    NC_000913   5088    .   A   T

This way, you can identify variants you expected to be called in your
samples, such as alleles targeted by MAGE oligonucleotides.
