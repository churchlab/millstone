import os

from conf.global_settings import *

try:
	from conf.local_settings import *
except ImportError:
	pass

# Use a default temp file root if none is set.
if 'TEMP_FILE_ROOT' not in locals():
    TEMP_FILE_ROOT = os.path.join(MEDIA_ROOT, 'tmp')

# Set the binaries here in case TOOLS_DIR is modified in local_settings.py
# TODO(gleb): What if users want to override specific binaries? Probably
# want a tools_settings.py.
SAMTOOLS_BINARY = '%s/samtools/samtools' % TOOLS_DIR
BGZIP_BINARY = '%s/tabix/bgzip' % TOOLS_DIR
FASTQC_BINARY = '%s/fastqc/fastqc' % TOOLS_DIR
VCFUNIQ_BINARY = '%s/freebayes/vcfuniq' % TOOLS_DIR
VCFSTREAMSORT_BINARY = '%s/freebayes/vcfstreamsort' % TOOLS_DIR

