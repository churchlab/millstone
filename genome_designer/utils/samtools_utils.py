"""
Wrappers for samtools scripts.

NOTE: Avoid adding any dependency on Django or Millstone to this module
as eventually we might want to extract it into it sown library.
The only exception to this for now is the link to the tool binary.
"""

import subprocess

from django.conf import settings


def run_mpileup(input_bam_file, input_ref_genome_fasta_path, output_path,
        coverage_only=True):
    with open(output_path, 'w') as fh:
        p_mpileup = subprocess.Popen([
                settings.SAMTOOLS_BINARY,
                'mpileup',
                '-f', input_ref_genome_fasta_path,
                input_bam_file
        ], stdout=subprocess.PIPE)

        if coverage_only:
            subprocess.check_call([
                    'cut',
                    '-f',
                    '-4'
            ], stdin=p_mpileup.stdout, stdout=fh)
        else:
            for line in p_mpileup.stdout:
                fh.write(line)
