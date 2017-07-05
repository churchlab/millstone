**************************
SNV Calling and Annotation
**************************

SNV calling is performed by Freebayes.

Annotation of variants is performed by SnpEff. The default arguments to SnpEff are specified in the code `here <https://github.com/churchlab/millstone/blob/972cf2e7c38d796ec49aebb77a1aec5742986f13/genome_designer/pipeline/variant_effects.py#L320>`_ and some can be overriden by modifying local_settings.py and restarting the Millstone web server and celery.
