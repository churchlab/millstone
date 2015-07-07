**DISCLAIMER: The current Millstone Amazon setup leaves your application
open to the web. Even though user accounts are password-protected,
certain uploaded and/or processed data is downloadable without
authentication if others "guess" the right urls. Realistically, this
shouldn't be a problem for most projects, but we're letting you know
just in case.**

This guide walks you through cloning the latest stable Amazon Machine
Image (AMI) configured with Millstone. The AMI will automatically set up
Millstone, all you need to do is clone it into an Amazon instance, start
the instance, and log in.

*All new users will want to use this guide. Docs for individuals wishing
to configure their instance or modify source code are coming soon.*

**Table of Contents**

-  `Before Reading this Guide <#before-reading-this-guide>`__
-  `Create an Amazon AWS Account <#create-an-amazon-aws-account>`__
-  `Cloning the AMI <#cloning-the-ami>`__
-  `Accessing your instance <#accessing-your-instance>`__
-  `In the browser <#in-the-browser>`__
-  `On the command line <#on-the-command-line-just-in-case>`__
-  `Using Millstone <#using-millstone>`__
-  `Registration <#registration>`__
-  `New Project and First
   Alignment <#new-project-and-first-alignment>`__
-  `Viewing Variants <#viewing-variants>`__
-  `Cast vs. Melted <#cast-vs-melted>`__
-  `Links <#links>`__
-  `Fields and Filtering <#fields-and-filtering>`__
-  `Examples <#examples>`__
-  `Viewing Variants <#viewing-variants>`__
-  `Cast vs. Melted <#cast-vs-melted>`__
-  `Links <#links>`__
-  `Fields and Filtering <#fields-and-filtering>`__
-  `Examples <#examples>`__
-  `Marginal Calls <#marginal-calls>`__
-  `Variant Sets <#variant-sets>`__
-  `Creating a blank set <#creating-a-blank-set>`__
-  `Uploading a set from a VCF
   file <#uploading-a-set-from-a-vcf-file>`__
-  `Troubleshooting <#troubleshooting>`__

Before Reading this Guide
-------------------------

**Please read the `Introduction to
Millstone <https://github.com/churchlab/churchlab.github.io/wiki/Introduction-to-Millstone>`__
before going through these steps.** That document will help you ensure
that Millstone is right for you, and that you have the right sort of
data.

Create an Amazon AWS Account
----------------------------

You need to create to an Amazon Web Services (AWS) account. `Brad
Chapman's getting started guide for
cloudbiolinux <https://github.com/chapmanb/cloudbiolinux/blob/master/doc/intro/gettingStarted_CloudBioLinux.pdf?raw=true>`__
has a solid first chapter with instructions on getting everything set
up.

Cloning the AMI
---------------

1. Login to https://console.aws.amazon.com/console/home and proceed to
   EC2. In the upper-right corner, be sure to select the N. Virgina
   region. We can't guarantee our AMI is visible outside of that region.
   From the EC2 dashboard, press ``Launch Instance``, which will take
   you into a Wizard to have you configure your instance.

2. In the Choose AMI tab, select Community AMIs in the left panel, then
   search for "millstone", or search directly using our beta release ami
   id: ami-c9746da0.

3. On the 'Choose instance type' tab, select an instance according to
   your needs. We recommend m3.medium (select General Purpose on the
   left). The number of vCPUs will determine how many genomes can be
   simultaneously aligned.

4. In 'Configure instance', the only setting we recommend changing is
   explicitly setting the Availability Zone (we always use
   ``us-east-1a``). You can only move EBS (Amazon hard drives) between
   instances in the same zone, so it'll make things easier to
   consistently make everything in the same zone.

5. In 'Add storage', increase the size of the root drive to the amount
   of space that you'll need. For bacterial genomes, about 2 GB per
   sample should be more than enough (i.e. 100 samples = 200 GB).

6. In 'Tag instance', fill in an informative value for the 'Name' key.
   We like the name to include the date it was created and a description
   of what the instance is running (e.g.
   ``2014_04_01_mutate_all_the_things``).

7. For security group, configure a group appropriate to your needs. Most
   users will want to create a security group with all of the following
   open. (*This will make your instance publicly visible to someone
   trying random EC2 IPs, but login is still required.*):

   -  All ICMP
   -  All TCP
   -  All UDP
   -  SSH

8. Continue to the final tab where you'll press 'Launch the instance'.
   Select or create a public/private key pair. If you create the key,
   download and save the private key, and put it somewhere safe (we
   suggest ``~/.ssh/``.) (*If you lose the private key there's no way to
   ssh back into your instance. You'll have to terminate it and create a
   new one.*)

It takes about 5-10 minutes for the instance to launch and all
bootstrapping to finish, after which your Millstone is ready to grind!

Accessing your instance
-----------------------

Go back to the `EC2 console Instances
page <https://console.aws.amazon.com/ec2/v2/home?#Instances:>`__ and
make sure you are in the correct region, using the dropdown in the top
right. The instance you created should be visible in the list. When it
is ready, its *Status Checks* column should say '2/2 checks passed'.

In the browser
~~~~~~~~~~~~~~

.. figure:: https://cloud.githubusercontent.com/assets/515076/6034315/591031d4-abef-11e4-87bd-d66286b31b15.png
   :alt: 

Select the instance from the list, and the info pane should appear below
the instance list. In the Description tab, the webpage URL to can be
found under **Public DNS**. The url should look like:
``ec2-xx-xx-xx-xx.compute-1.amazonaws.com``

It may take some time for your instance to initialize. Wait until all
status checks are completed before attempting to log in. If the server
doesn't come up, it might still be loading.

On the command line (just in case)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It should not be necessary at the moment, but if you need to SSH into
the server, the command is:

::

    ssh -i ~/.ssh/your-key.pem ubuntu@ec2-xx-xx-xx-xx.compute-1.amazonaws.com

(This assumes you put the private key you generated in ``~/.ssh/``). If
permissions fail on your key, ``chmod`` the key's permissions to 700.

Using Millstone
---------------

Registration
~~~~~~~~~~~~

Once Millstone is installed, you should be greeted with the Millstone
logo and a login/register page. Register a user with a login, email, and
password. *Currently, we only allow one user per instance. After the
first user is registered, registration is closed.* **Don't forget your
username and password, as there is currently no 'reminder'
functionality.** (The only way to change your password at present is to
do so through the Django shell, using methods available on the Django
auth model ``User``.)

New Project and First Alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
^^^^^^^^^^^^^^^^

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
^^^^^^^

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
filenames (no path required). Here is an example:

::

    Sample_Name Read_1_Filename Read_2_Filename
    sample01    sample01_fwd.fq.gz  sample01_rev.fq.gz
    sample02    sample02_fwd.fq.gz  sample02_rev.fq.gz

*NOTE: Millstone can work with ``gzip``-ed FASTQ files, and they will be
faster to upload.*

Once you upload the template, it will list the samples awaiting upload:

.. figure:: https://cloud.githubusercontent.com/assets/515076/6034627/fb35f306-abf2-11e4-947c-e14664e9c804.png
   :alt: 

You can then upload the individual files matching the filenames in the
template.

Alignment Settings
^^^^^^^^^^^^^^^^^^

By default, Millstone treats all samples as diploid. This allows
ambiguous variants to be called as heterozygous. You can choose to keep
all of these ambiguous variants, to keep only those where at least some
samples are called as non-ambiguous, or throw away ambiguous variants
all together. If you have many samples, we suggest the latter two
options to keep the database size manageable.

Submit Alignment
^^^^^^^^^^^^^^^^

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

Viewing Variants
----------------

All sample, reference genome, and alignments are listed in the *Data*
view (the toggle switch in the top left). Clicking over to the *Analyze*
view will allow you to filter through multi-sample variants and view
their aligned reads. Use the dropdowns on the left below the
Data/Analyze toggle to select your alignment and your reference genome
and choose *Variants*. Once the alignments are complete, you should see
a list of all variants that have been identified across all samples.

Cast vs. Melted
~~~~~~~~~~~~~~~

There are two ways to view variants.

**Cast**: Cast displays a summary row for one variant across all
samples. You can see how many samples the variant is present in, as well
as the variant's effects.

**Melted**: 'Melting' the view shows one row for every combination of
sample and and variant. It essentially multiplies the rows by the number
of samples, so you can see data specific to individual samples. If a
variant is not called in a sample, it's *Alt* column will be blank.

Links
~~~~~

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
~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~

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
~~~~~~~~~~~~

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

Troubleshooting:
----------------

-  I can log in via SSH but the web interface doesn't load!

    You've probably forgotten to allow access to your instance through
    web interfaces. This can be fixed by adding the following
    connections to your security group: \* All ICMP \* All TCP \* All
    UDP You can do this by going to the Network & Security -> Security
    Groups section of the EC2 dashboard and editing the security group
    that you created in your instance. If you've forgotten this can be
    found in the main instance dash on the far right under security
    groups. Click on that and you should be able to edit inbound rules
    by right clicking on the Group ID

-  I've managed to load the webpage but get a 502 bad gateway error!

    Millstone is probably loading up, try again in a few minutes.

-  Registration is closed.

    Only one user is allowed to register (as soon as the server boots
    up), and afterwards registration is closed.

-  Millstone just sits there after importing a template file.

    This could be any number of things. If your template file is
    formatted correctly, it could be a completely out of space error, so
    check that you've got room on your drive containing Millstone. File
    formatting is often the biggest problem in this stage, so be careful
    that you've escaped spaces in file names.

-  I want to make sure everything's going right, where can I find the
   logs?

    The logs are by default at /var/log/supervisor
