"""
Methods related to importing data.
"""

import os
import shutil
import re

from main.models import Dataset
from main.models import ReferenceGenome
from settings import PWD

from Bio import SeqIO


IMPORT_FORMAT_TO_DATASET_TYPE = {
    'fasta': Dataset.TYPE.REFERENCE_GENOME_FASTA,
    'genbank': Dataset.TYPE.REFERENCE_GENOME_GENBANK
}

def import_reference_genome_from_local_file(project, label, file_location,
        import_format):
    """Creates a ReferenceGenome associated with the given Project.

    Args:
        project: The Project we're storing everyting relative to.
        label: The human-readable label for the ReferenceGenome.
        file_location: Location of the genome on the server.
        import_format: Must be 'fasta' or 'genbank'.
    """
    # Validate the input.
    assert import_format in ['fasta', 'genbank']

    # Validate the file.
    assert os.path.exists(file_location), "File doesn't exist."

    # Validate the input by parsing it with BioPython, while also
    # counting the number of chromosomes.
    num_chromosomes = 0
    num_bases = 0
    for genome_record in SeqIO.parse(file_location, import_format):
        num_chromosomes += 1
        num_bases += len(genome_record)

    # Create the ReferenceGenome object.
    reference_genome = ReferenceGenome.objects.create(
            project=project,
            label=label,
            num_chromosomes=num_chromosomes,
            num_bases=num_bases)

    # Copy the source file to the ReferenceGenome data location.
    dataset_type = IMPORT_FORMAT_TO_DATASET_TYPE[import_format]
    _copy_and_add_dataset_source(reference_genome, dataset_type,
            dataset_type, file_location)
    return reference_genome

def import_samples_from_targets_file(project, targets_file):
    """Uses the uploaded targets file to add a set of samples to the project.
    We need to check each line of the targets file for consistency before we
    do anything, however. 
    
    It writes a copy of the uploaded targets file to a temporary file
    
    Args:
        project: The project we're storing everything relative to
        targets_file: The UploadedFile django object that holds the targets
            in .tsv format. 
    """
    
    # The targets file shouldn't be over 1 Mb ( that would be ~3,000 genomes)
    assert targets_file.size < 1000000, "Targets file is too large."
    
    #Check the header first, and make sure it has the required columns in 
    # the right order.
    target_header = targets_file.readline().rstrip().split('\t')
    
    assert len(target_header) >= 6, "Bad header. Were columns removed?"
    
    required_header = ['Sample_Name','Plate_or_Group','Well','Read_1_Path','Read_2_Path','Parent_Samples']
    for col, check in zip(target_header[0:len(required_header)], required_header):
        assert col == check, "Header column '%s' is missing or out of order." 
    
    optional_header = target_header[len(required_header):]
    field_names = required_header + optional_header 
    
    #Now go through every row and make a dict. 
    #All paths should be correct and all values should be alphanumeric.
    #Also, we need to make sure that we have access to that path
    rows = []
    row_count = -1
    for row in targets_file:
        row_dict = {}
        row_count += 1
        #skip the first row, which is the header
        if row_count == 0: continue
                
        fields = row.rstrip().split('\t')
        
        assert len(fields) == len(field_names), "Row %d has the wrong length." % row_count
        
        for field_name, field in zip(field_names, fields): 
            
            print "%d : %s : %s" % (row_count, field_name, field)
            
            if 'Path' not in field_name:
                #make sure each field is alphanumeric only
                assert re.match('^[\. \w-]*$', field) is not None, 'Only alphanumeric characters and spaces are allowed, except for the paths.\n(Row %d, "%s")' % (row_count, field_name)
            else:
                #if it is a path, then try to open the file and read one byte
                #replace the string '$GD_ROOT with the project path, so we can use the test data
                field = field.replace('$GD_ROOT', PWD)
                
                with open(field, 'rb') as test_file:
                    try:
                        test_file.read(8)
                    except:
                        raise
            
            row_dict[field_name] = field
            
        rows.append(row_dict)
    
    #TODO: use this list of dicts to update the database, create new sample dirs,
    #      and copy the data to a location within GD.     
    for row in rows:
        print row

def _copy_and_add_dataset_source(entity, dataset_label, dataset_type,
        original_source_location):
    """Copies the dataset to the entity location and then adds as
    Dataset.

    Returns the Dataset object.
    """
    source_name = os.path.split(original_source_location)[1]
    dest = os.path.join(entity.get_model_data_dir(), source_name)
    if not original_source_location == dest:
        shutil.copy(original_source_location, dest)
    return _add_dataset_to_genome(
            entity, dataset_label, dataset_type, dest)


def _add_dataset_to_genome(genome, dataset_label, dataset_type,
        filesystem_location):
    """Helper function for adding a Dataset to a ReferenceGenome.

    Returns the Dataset object.
    """
    dataset = Dataset.objects.create(
            label=dataset_label,
            type=dataset_type,
            filesystem_location=filesystem_location)
    genome.dataset_set.add(dataset)
    genome.save()
    return dataset
