from conf.global_settings import *
from conf.local_settings import *

TEMP_FILE_ROOT = os.path.join(MEDIA_ROOT, 'tmp')


# Set the binaries here in case TOOLS_DIR is modified in local_settings.py
# TODO(gleb): What if users want to override specific binaries? Probably
# want a tools_settings.py.
SAMTOOLS_BINARY = '%s/samtools/samtools' % TOOLS_DIR
BGZIP_BINARY = '%s/tabix/bgzip' % TOOLS_DIR
