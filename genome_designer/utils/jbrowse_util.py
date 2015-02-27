"""
Utility methods for creating JBrowse config files to allow data
to be viewed using JBrowse.

## Note about JSON Track Concurrency Fix:

We encountered the possibility of multiple workers wanting to
add/remove/change tracks simultaneously in the tracksList.json file. This
results in race conditions where one worker overwrites the changes of another,
causing missing tracks.

To fix this problem, we write to a directory structure inside the jbrowse
directory called ./indiv_tracks, where tracks are stored in subdirectories
under the ref genome by their label.

`jbrowse/indiv_tracks/{track['label']}

When an alignment completes, compile_tracklist_json() is called, which gathers
up all the tracks for a reference genome and updates the 'tracks' field of the
json tracklist. It also needs to symlink various track subdirectories into the
main directory. It assumes that there won't be anything to overwrite (like two
{trackname}/seq dirs).

We perform this compile_tracklist_json() every time a link is asked for.

"""

from distutils.dir_util import mkpath
import glob
import json
import os
import shutil
import subprocess

from main.model_utils import get_dataset_with_type
from main.models import Dataset
from utils import merge_nested_dictionaries
from settings import JBROWSE_BIN_PATH
from settings import JBROWSE_DATA_URL_ROOT
from settings import JBROWSE_GBK_TYPES_TO_DISPLAY
from settings import S3_BUCKET
from settings import TOOLS_DIR

TABIX_BINARY = '%s/tabix/tabix' % TOOLS_DIR

# TODO: Figure out better place to put this.
# JBrowse requires the symlink path to exist. See settings.py
# comments for more info.

# assert os.path.exists(JBROWSE_DATA_SYMLINK_PATH), (
#         "%s does not exists. You may need to symlink it." %
#                 JBROWSE_DATA_SYMLINK_PATH)

def get_tracklist_json(reference_genome, concurrent_id=None):
    """
    Get trackList for jbrowse
    """

    jbrowse_path = reference_genome.get_jbrowse_directory_path()

    # If this is a concurrent read, read from individual tracks
    # subdir, named by a concurrent_id. Update is OK if dir/id is
    # there.
    if concurrent_id:
        jbrowse_path = os.path.join(
                jbrowse_path, 'indiv_tracks', concurrent_id)
        assert os.path.exists(jbrowse_path), (
                'concurrent_id for json track does not exist: %s' % (
                    jbrowse_path))

    # Load the old tracks
    json_track_fn = os.path.join(jbrowse_path, 'trackList.json')
    return json.loads(open(json_track_fn).read())

def write_tracklist_json(reference_genome, dictionary, concurrent_id=None):
    """
    Overwrite trackList for jbrowse

    """
    jbrowse_path = reference_genome.get_jbrowse_directory_path()

    # If this is a concurrent write, write to an individual tracks
    # subdir, named by a concurrent_id. Update is OK if dir/id is
    # there.
    if concurrent_id:
        jbrowse_path = os.path.join(
                jbrowse_path, 'indiv_tracks', concurrent_id)

        # distutils makes intermediate directories only if necessary,
        # os.makedirs isn't as forgiving
        mkpath(jbrowse_path)

    json_track_fn = os.path.join(jbrowse_path, 'trackList.json')

    # Create a new or overwrite an old tracklist
    with open(json_track_fn, 'w') as json_track_fh:
        json_track_fh.write(json.dumps(dictionary))

def compile_tracklist_json(reference_genome):
    """
    Gathers all the individual tracks in the ./indiv_tracks
    directory and creates a new 'tracks' listing, keeping
    each track label unique.
    """
    jbrowse_path = reference_genome.get_jbrowse_directory_path()
    indiv_tracks_path = os.path.join(jbrowse_path,'indiv_tracks')
    track_files = glob.glob(os.path.join(indiv_tracks_path,'*','trackList.json'))

    # a dictionary of tracks by label. We assume here that all
    # tracks are unique by 'label' (which should really be called
    # 'key', since its a machine-readable field)
    track_dict = {}
    consolidated_track_list = {}

    for track_fn in track_files:

        this_track_list = json.loads(open(track_fn).read())

        # First, pop out off 'tracks' list from the json, and
        # add them to the track_dict by label
        for track in this_track_list.pop('tracks', []):
            track_key = track['label']

            track_dict[track_key] = track

            # prepend the indiv_track subdir to any relative paths
            if 'urlTemplate' in track:
                if track['urlTemplate'].startswith('/'): continue

                track['urlTemplate'] = os.path.join(
                    'indiv_tracks',track_key, track['urlTemplate'])


        # Finally, symlink any subdirs (seq, etc) into the root.
        track_dir = os.path.dirname(track_fn)

        for subdir in os.listdir(track_dir):
            abs_subdir = os.path.join(track_dir, subdir)
            if not os.path.isdir(abs_subdir): continue

            try:
                shutil.copytree(
                    src= abs_subdir,
                    dst= os.path.join(jbrowse_path,subdir),
                    symlinks=True
                )
            except Exception as e:
                print e
                print 'Skipping dir {:s} because it already exists.'.format(
                        subdir)

        # (We're assuming here that any overwriting that individual
        # files do of these fields is not important. The only field
        # that I am aware of is 'formatVersion', and is 1 for all.)
        consolidated_track_list = merge_nested_dictionaries(
                consolidated_track_list, this_track_list)

    # Add back all the tracks
    consolidated_track_list['tracks'] = track_dict.values()

    # touch the other tracklist format, make it empty
    open(os.path.join(jbrowse_path, 'tracks.conf'), 'a',).close()

    write_tracklist_json(reference_genome, consolidated_track_list)


def prepare_jbrowse_ref_sequence(reference_genome, **kwargs):
    """Prepare the reference sequence and place it in the ref_genome dir.

    This implicitly creates the config directory structure for this reference
    genome. Tracks added in the future are added relative to this reference
    genome.

    The implemenation of this method is a light wrapper around
    jbrowse/bin/prepare-refseqs.pl.
    """
    PREPARE_REFSEQS_BIN = os.path.join(JBROWSE_BIN_PATH, 'prepare-refseqs.pl')

    # First ensure that the reference genome exists. If it fails, try to
    # convert from genbank, then give up.
    reference_fasta = get_dataset_with_type(
            reference_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_FASTA).get_absolute_location()

    # Next, ensure that the jbrowse directory exists.
    reference_genome.ensure_jbrowse_dir()
    jbrowse_path = os.path.join(
            reference_genome.get_jbrowse_directory_path(),
            'indiv_tracks',
            'DNA')

    # Now run prepare-refseqs.pl to get the ReferenceGenome in.
    subprocess.check_call([
        PREPARE_REFSEQS_BIN,
        '--fasta', reference_fasta,
        '--out', jbrowse_path,
    ])

    json_tracks = get_tracklist_json(reference_genome, 'DNA')

    # DNA track should be the first track
    dna_track = json_tracks['tracks'][0]
    assert dna_track['type'] == 'SequenceTrack'

    # Get rid of translation and reverse strand
    dna_track.update({
        "showForwardStrand": True,
        "showReverseStrand": False,
        "showTranslation": False
        })

    write_tracklist_json(reference_genome, json_tracks, 'DNA')

def add_genbank_file_track(reference_genome, **kwargs):
    """
    Jbrowse has the ability to make tracks out of genbank files. This
    takes the genbank file from a reference_genome object and attempts to
    make such a track and then add it to the track list.
    """
    FLATFILE_TRACK_BIN = os.path.join(JBROWSE_BIN_PATH, 'flatfile-to-json.pl')

    reference_gbk = get_dataset_with_type(
            reference_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_GENBANK).get_absolute_location()

    jbrowse_path = reference_genome.get_jbrowse_directory_path()

    reference_gff = get_dataset_with_type(
            reference_genome,
            type=Dataset.TYPE.REFERENCE_GENOME_GFF).get_absolute_location()

    json_update_fields = {
        'style': {
            'label': 'name,CDS,gene',
            'description': 'note,function,gene_synonym',
            'color': '#43a110'
        }
    }

    genbank_json_command = [
        FLATFILE_TRACK_BIN,
        '--gff', reference_gff,
        '--out', os.path.join(jbrowse_path,'indiv_tracks','gbk'),
        '--type', JBROWSE_GBK_TYPES_TO_DISPLAY,
        '--autocomplete','all',
        '--trackLabel','gbk',
        '--key',"Genome Features",
        '--trackType',"CanvasFeatures",
        #'--getSubfeatures',
        #'--className','transcript',
        #'--subfeatureClasses', "{\"CDS\":\"transcript-CDS\"}"
    ]

    subprocess.check_call(genbank_json_command)

    # Finally, manually update tracklist json with style info
    tracklist_json = get_tracklist_json(reference_genome, 'gbk')

    for i, track in enumerate(tracklist_json['tracks']):
        if track['key'] == 'Genome Features':
            tracklist_json['tracks'][i] = merge_nested_dictionaries(
                    track, json_update_fields)

    write_tracklist_json(reference_genome, tracklist_json, 'gbk')


def add_vcf_track(reference_genome, alignment_group, vcf_dataset_type):
    """Adds a vcf track to JBrowse for this vcf.

    See JBrowse Docs:
        http://gmod.org/wiki/JBrowse_Configuration_Guide#Example_VCF-based_Variant_Track_Configuration
    """
    # Get the vcf file location from the the Dataset of the genome
    # keyed by the alignment_type.
    vcf_dataset = get_dataset_with_type(alignment_group,
            vcf_dataset_type)

    vcf_dataset = _vcf_to_vcftabix(vcf_dataset)

    if reference_genome.project.is_s3_backed():
        urlTemplate = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
            vcf_dataset.filesystem_location.strip("/jbrowse"))
        urlTemplate_idx = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
            vcf_dataset.filesystem_idx_location.strip("/jbrowse"))
    else:
        urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
            vcf_dataset.filesystem_location)
        urlTemplate_idx = os.path.join(JBROWSE_DATA_URL_ROOT,
            vcf_dataset.filesystem_idx_location)

    label = vcf_dataset.internal_string(alignment_group)
    key = "{:s} {:s} SNVs".format(vcf_dataset.type,alignment_group.label)

    # Build the JSON object.
    raw_dict_obj = {
        'tracks' : [{
            "label"         : label,
            "key"           : key,
            "storeClass"    : "JBrowse/Store/SeqFeature/VCFTabix",
            "urlTemplate"   : urlTemplate,
            "tbiUrlTemplate": urlTemplate_idx,
            'category'      : 'VCF Tracks',
            "type"          : "JBrowse/View/Track/HTMLVariants"
        }]
    }

    write_tracklist_json(reference_genome, raw_dict_obj, label)


def _vcf_to_vcftabix(vcf_dataset):
    """Compresses and indexes a vcf using samtools tabix.

    Creates a new Dataset model instance for this compressed version, with the
    same related objects (e.g. pointing to the same AlignmentGroup). The
    Dataset is flagged as compressed, indexed, etc.

    Args:
        vcf_dataset: Dataset pointing to a vcf, or its compressed version.
            Index may or may not exist.

    Returns:
        Dataset that points to compressed version of input vcf_dataset, if it
        wasn't compressed already. The index file is asserted to exist for this
        compressed Dataset.
    """
    ### This function has two steps:
    ###     1. Get or create compressed Dataset.
    ###     2. Create index if it doesn't exist.

    ### 1. Get or create compressed Dataset.
    if vcf_dataset.is_compressed():
        compressed_dataset = vcf_dataset
    else:
        # Check for existing compressed version using related model.
        # Assume that the first model will do.
        related_model = vcf_dataset.get_related_model_set().all()[0]
        compressed_dataset = get_dataset_with_type(
                entity=related_model,
                type=vcf_dataset.type,
                compressed=True)
        # If there is no compressed dataset, then make it
        if compressed_dataset is None:
            compressed_dataset = vcf_dataset.make_compressed('.bgz')

    ### 2. Create index if it doesn't exist.
    if compressed_dataset.filesystem_idx_location == '':

        # Set the tabix index location
        compressed_dataset.filesystem_idx_location = (
                compressed_dataset.filesystem_location + '.tbi')
        compressed_dataset.save()

        # Make tabix index
        subprocess.check_call([
            TABIX_BINARY, '-f',
            '-p', 'vcf',
            compressed_dataset.get_absolute_location()
        ])

    # Make sure the index exists, whether created now or previously.
    assert compressed_dataset.filesystem_idx_location == (
            compressed_dataset.filesystem_location + '.tbi'), (
            'Tabix index file location is not correct.')
    assert os.path.exists(
            compressed_dataset.get_absolute_idx_location()), (
            'Tabix index file does not exist on filesystem.')

    return compressed_dataset


def add_bed_file_track(reference_genome, sample_alignment, dataset):
    """ Add a bed file to Jbrowse, like that created for CallableLoci.
        Pass in the dataset model object directly.
    """

    FLATFILE_TRACK_BIN = os.path.join(JBROWSE_BIN_PATH, 'flatfile-to-json.pl')

    bed_dataset_location = dataset.get_absolute_location()

    jbrowse_path = reference_genome.get_jbrowse_directory_path()

    # doing label as ES_AG because SA isn't currently used in the variant view
    label = dataset.internal_string(sample_alignment.experiment_sample) + '_' + (
            str(sample_alignment.alignment_group.uid))
    key = ':'.join([
            sample_alignment.alignment_group.label,
            dataset.external_string(sample_alignment.experiment_sample)])

    bed_json_command = [
        FLATFILE_TRACK_BIN,
        '--bed', bed_dataset_location,
        '--out', os.path.join(jbrowse_path,'indiv_tracks',label),
        '--trackLabel',label,
        '--key', key,
        '--trackType',"CanvasFeatures"
    ]

    subprocess.check_call(bed_json_command)

    # Finally, manually update category for json
    tracklist_json = get_tracklist_json(reference_genome, label)

    for i, track in enumerate(tracklist_json['tracks']):
        track['category'] = 'BED Features'

    write_tracklist_json(reference_genome, tracklist_json, label)


def add_bam_file_track(reference_genome, sample_alignment, alignment_type):
    """Update the JBrowse track config file, trackList.json, for this
    ReferenceGenome with a track for the given sample_alignment and alignment_type.
    """
    # Get the bam file location from the the Dataset of the genome
    # keyed by the alignment_type.
    bam_dataset = get_dataset_with_type(sample_alignment, alignment_type)

    # Figure out the url that JBrowse would use to show the data, e.g.:
    #     /jbrowse/gd_data/projects/58a62c7d/genomes/8dc829ec/align.bam
    # urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
    #         bam_dataset.filesystem_location)

    # NOTE: We should construct bam file urls using project.get_client_jbrowse_link()
    #       rather than checking S3 flag here.
    if reference_genome.project.is_s3_backed():
        urlTemplate = os.path.join('http://%s.s3.amazonaws.com/' % S3_BUCKET,
            bam_dataset.filesystem_location.strip("/jbrowse"))
    else:
        urlTemplate = os.path.join(JBROWSE_DATA_URL_ROOT,
            bam_dataset.filesystem_location)

    # doing label as ES_AG because SA isn't currently used in the variant view
    label = bam_dataset.internal_string(sample_alignment.experiment_sample) + '_' + (
            str(sample_alignment.alignment_group.uid))
    key = ':'.join([
            sample_alignment.alignment_group.label,
            bam_dataset.external_string(sample_alignment.experiment_sample)])

    # Build the JSON object.
    raw_dict_obj = {
        'tracks' : [
        {
            'storeClass': 'JBrowse/Store/SeqFeature/BAM',
            'urlTemplate': urlTemplate,
            'label': label,
            'type': 'JBrowse/View/Track/Alignments2',
            'chuckSizeLimit': 10000000, #double the default chunk size
            'key': key,
            'category': 'BAM Tracks',
            'style' : {
                'className': 'alignment',
                'arrowheadClass': 'arrowhead',
                'labelScale': 100
            }
        }
    ]}
    write_tracklist_json(reference_genome, raw_dict_obj, label)

    # Also add a snp coverage track.
    snp_coverage_label = bam_dataset.internal_string(
            sample_alignment.experiment_sample) + '_COVERAGE_' + (
            str(sample_alignment.alignment_group.uid))

    snp_coverage_key = key + ' Coverage'
    coverage_raw_dict_obj = {
        'tracks' : [
        {
            'storeClass': 'JBrowse/Store/SeqFeature/BAM',
            'urlTemplate': urlTemplate,
            'label': snp_coverage_label,
            'type': 'JBrowse/View/Track/SNPCoverage',
            'category': 'Coverage Tracks',
            'key': snp_coverage_key
        }
    ]}
    write_tracklist_json(reference_genome, coverage_raw_dict_obj,
            snp_coverage_label)

def get_tracks_for_entity(reference_genome,
        alignment_group=None,
        sample_alignment=None):
    """
    Get a list of jbrowse track labels (i.e. the internal track names)
    that we want to show. Stored in a tuple, with the first string as the
    category.
    """
    raise NotImplementedError
    # This might be a good idea later.

    # track_list = defaultdict(dict)

    # # Reference Genome

    # # always have DNA and gbk if available
    # track_list['DNA']('DNA')
    # if reference_genome.is_annotated:
    #     track_list.append('gbk')

    # # Do all alignment groups if None, else get just the one:
    # if alignment_group is not None:
    #     alignment_groups = reference_genome.alignmentgroup_set.all():
    # else:
    #     alignment_groups = [alignment_group]

    # for alignment_group in alignment_groups:

    #     #VCF
    #     track_list.append('vcf',
    #                 '_'.join([
    #         str(alignment_group.uid),
    #         Dataset.TYPE.VCF_FREEBAYES_SNPEFF]))




