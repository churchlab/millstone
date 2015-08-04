********************
Setting up Millstone
********************

Using Millstone on Amazon AWS
=============================

**Using Millstone via AWS is the preferred option for most users.** We have pre-configured
a Millstone installation into an Amazon Machine Image (AMI). This means you can sidestep
all of the dependency installation, configuration, etc.

**DISCLAIMER: The current Millstone Amazon setup leaves your application
open to the web. Even though user accounts are password-protected,
certain uploaded and/or processed data is downloadable without
authentication if others "guess" the right urls. Realistically, this
shouldn't be a problem for most projects, but we're letting you know
just in case.**

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
   id: ami-4b8d5220.

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

Using Millstone locally
=======================

It is also possible to use Millstone locally on Mac OSX and Linux. **Local installation is meant for advanced users only**. It requires the manual installation and configuration of various dependencies, and requires root access. You can our find local installation guide in the `Millstone github readme <https://github.com/churchlab/millstone#installation>`__.
