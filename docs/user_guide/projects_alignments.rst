******************************
Projects and Alignments
******************************

Registering a new user
======================

Once Millstone is installed, you should be greeted with the Millstone
logo and a login/register page. Register a user with a login, email, and
password. *Currently, we only allow one user per instance. After the
first user is registered, registration is closed.* **Don't forget your
username and password, as there is currently no 'reminder'
functionality.** (The only way to change your password at present is to
do so through the Django shell, using methods available on the Django
auth model ``User``.)

Creating a new project
======================

Once you register, you can create a new project, and you will the
prompted to give it a short name. Afterwards, you will be taken to the
create alignment screen. There are 5 steps, each with a tab in the top
bar. Choose a name for your first alignment, which will pair a reference
genome with a set of samples to align. One project can have multiple
alignments.

*Note: If you have many/large samples, and would prefer to upload files
via the command line instead of the browser, see `this
guide <https://github.com/churchlab/churchlab.github.io/wiki/Manual-Data-Upload-to-Server>`__.*

Reference Genome
----------------

Select the Reference Genome tab, and click the green 'New' button. You
can select a reference genome from NCBI or upload a custom reference.

*Note: If you use a FASTA there will be no variant annotation
information, so Genbank is recommended if you have one.*

**Load file from NCBI**: Simply fill in the accession number (for
instance `U00096.2 <http://www.ncbi.nlm.nih.gov/nuccore/U00096.2>`__ for
E. coli) and give the reference genome a name. If you'd like to use a
custom reference genome, you can upload a file from your desktop. You
can check to make sure you've got the right accession number by
comparing your genome's size to the number of nucleotides present in the
reference genome.

**Upload through browser**: If you have a local file with your genome,
you can upload it with this option. If you have a large cassette
insertion or plasmid you would also like to align, you can edit the
FASTA/Genbank file to insert it into the genome using a tool like
Benchling or Geneious (in the case of a cassette insertion), or add it
as a separate chromosome (an additional FASTA or GenBank record in the
same file).

Finally, select the checkbox next to the uploaded genome to mark it as
your reference.

Samples
-------

Once that's done, move on to the samples tab. Each genome sample you
upload must consist of a pair of forward and reverse FASTQ files. You
can either upload samples through the browser, or you can upload them in
batch to the server using a the command line via ``scp``. The command
line approach is better for large numbers of samples, but is more
complicated. It is detailed in the *Manual Upload* section at the bottom
of this guide.

Open the upload samples dialog via the green 'New' button, then choose
'Batch Upload through browser...'. In order to upload samples through
the browser, you must first register samples to be uploaded by filling
out a spreadsheet template with sample labels and corresponding data
filenames (no path required). *Fields must be separated by tabs.* Here is an example::

    Sample_Name Read_1_Filename Read_2_Filename
    sample01    sample01_fwd.fq.gz  sample01_rev.fq.gz
    sample02    sample02_fwd.fq.gz  sample02_rev.fq.gz

*NOTE: Millstone can work with ``gzip``-ed FASTQ files, and they will be
faster to upload.*

You can also include additional columns as per-sample metadata, like growth
rates, plate and well, strain parentage, etc. Here is an example:

.. literalinclude:: ../_examples/targets_example_with_metadata.tsv

Once you upload the template, it will list the samples awaiting upload:

.. figure:: https://cloud.githubusercontent.com/assets/515076/6034627/fb35f306-abf2-11e4-947c-e14664e9c804.png
   :alt:

You can then upload the individual files matching the filenames in the
template.

Alignment Settings
------------------

By default, Millstone treats all samples as diploid. This allows
ambiguous variants to be called as heterozygous. You can choose to keep
all of these ambiguous variants, to keep only those where at least some
samples are called as non-ambiguous, or throw away ambiguous variants
all together. If you have many samples, we suggest the latter two
options to keep the database size manageable.

Submit Alignment
----------------

Finally! Click the *Run Alignment* button in the last tab to start the
alignment. Depending on your genome size, number of samples, and the
size of the instance you chose, this could take time. You can see how
individual sample alignments are progressing by clicking on the name of
the alignment in the label column of the Alignments view. Every sample
will have an *output log* link and a Job Status.

After the individual samples are done aligning, the Alignment status
will change to ``VARIANT_CALLING`` as variants across all samples are
called in aggregate. Once this step has completed, then the Alignment
status will read ``COMPLETED`` and you can switch to the *Analyze* view
to examine the called variants.