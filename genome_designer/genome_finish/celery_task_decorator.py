from functools import wraps
import os
import traceback

from main.models import ExperimentSampleToAlignment


def set_assembly_status(sample_alignment, status, force=False):

    # Make sure assembly status is not FAILED
    if not force:
        assert sample_alignment.data['assembly_status'] != (
                ExperimentSampleToAlignment.ASSEMBLY_STATUS.FAILED)

    # Set assembly status for UI
    sample_alignment.data['assembly_status'] = (
                ExperimentSampleToAlignment.ASSEMBLY_STATUS.ASSEMBLING)
    sample_alignment.save()


def report_failure_stats(file_name):
    """Decorator that writes to file the traceback and exception of the
    decorated function, which must have one argument that is an instance
    of the ExperimentSampleToAlignment model class that it will set the
    assembly status of to failed
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as exc:

                sample_alignment_args = []
                for arg in args:
                    if isinstance(arg, ExperimentSampleToAlignment):
                        sample_alignment_args.append(arg)

                assert len(sample_alignment_args) == 1
                sample_alignment = sample_alignment_args[0]

                # Set assembly status to FAILED
                sample_alignment.data['assembly_status'] = (
                        ExperimentSampleToAlignment.ASSEMBLY_STATUS.FAILED)
                sample_alignment.save()

                # Write exception with traceback to file
                tb = traceback.format_exc()
                file_path = os.path.join(
                        sample_alignment.get_model_data_dir(),
                        file_name)
                with open(file_path, 'w') as fh:
                    fh.write('tracback:%s\nexception:%r' % (tb, exc))

                # Raise the excepted exception
                raise exc

        return wrapper
    return decorator
