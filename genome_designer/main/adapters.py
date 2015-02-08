"""Adapter - Converts django models to Front-end format.

    Calls on the model's get_field_order method to know which fields that
    are available for display and what order. Returns a config json
    that dumps the field config and the json data for model/table.
"""

import datetime
import json
import string

from django.db.models import Model
from django.db.models.query import QuerySet

from main.constants import UNDEFINED_STRING
from main.models import ExperimentSample

OBJ_LIST = 'obj_list'

# Fields that are not explicitly displayed but are used on the frontend
# to update views of particular columns.
# TODO: Make each field encapsulate its own special display logic.
HIDDEN_PAIR_FIELDS = [
    'href',
    'uid',
]


def adapt_model_to_frontend(model, filters={}, obj_list=None, **kwargs):
    """Converts django models to frontend format.

    Calls on the model's get_field_order method to know which fields that
    are available for display and what order.

    Args:
        model: The django model that we are adapting
        filters: The filter conditions on the model (i.e. select
            which instances of the model to select.)
        obj_list: If provided, use this as the set of objects to filter.
            NOTE: I actually want to permanently change the interface to this
            rather than doing filtering here.

    Returns:
        A config json that dumps the field config and the json data for
        model/table display.
    """

    # Fetch all objects in this model and package as json to be rendered into
    # the dom. The data will be displayed to the user, e.g. via the javascript
    # DataTables component.

    # Get all objects that pass the filter.
    if obj_list is None:
        obj_list = model.objects.filter(**filters)

    # A list of dicts with object data, where each dict is one object
    # and all the fields required for front-end display.
    fe_obj_list = [
            adapt_model_instance_to_frontend(obj, **kwargs)
            for obj in obj_list]

    # Get a list of fields required for displaying the objects, in the order
    # in which they should be displayed.
    field_dict_list = model.get_field_order(**kwargs)

    # Each field is a dict with two keys, 'field' for field name and 'verbose'
    # for display name. Get each. If 'verbose' is missing, then make verbose
    # be the field with _'s turned to spaces and Title Cased.
    field_list = [fdict['field'] for fdict in field_dict_list]

    # Get the verbose field names, which will be used as column headers.
    def _get_verbose(fdict):
        if 'verbose' in fdict:
            return fdict['verbose']
        else:
            return string.capwords(fdict['field'],'_').replace('_',' ')
    field_verbose_names = [_get_verbose(fdict) for fdict in field_dict_list]

    # A list of dicts containing the order of each column and the field titles
    # for each column, used for configuring jquery.datatables.js
    obj_field_config = [{
        'mData': name,
        'sTitle': verbose_name
    } for (name, verbose_name) in zip(field_list, field_verbose_names)]

    # Package the result.
    return json.dumps({
        OBJ_LIST: fe_obj_list,
        'field_config': obj_field_config
    })


def adapt_experiment_samples_to_frontend(filters={}, obj_list=None, **kwargs):
    """ The sample metadata fields require their own custom adapter. """
    # Get all objects that pass the filter.
    if obj_list is None:
        obj_list = ExperimentSample.objects.filter(**filters).order_by('label')

    json_fields = {}
    for obj in obj_list:
        json_field_dicts = dict(
                [(key,{'field':key}) for key in obj.data.keys()])
        json_fields.update(json_field_dicts)

    # A list of dicts with object data, where each dict is one object
    # and all the fields required for front-end display.
    fe_obj_list = []
    for obj in obj_list:
        # default to empty string
        obj_json_fields = dict((field, '') for field in json_fields)
        obj_json_fields.update(obj.data)
        fe_obj_list.append(adapt_model_instance_to_frontend(obj,
                field_info= obj_json_fields,
                **kwargs))


    # Get a list of fields required for displaying the objects, in the order
    # in which they should be displayed.
    field_dict_list = ExperimentSample.get_field_order(**kwargs)
    field_dict_list.extend(json_fields.values())

    # Each field is a dict with two keys, 'field' for field name and 'verbose'
    # for display name. Get each. If 'verbose' is missing, then make verbose
    # be the field with _'s turned to spaces and Title Cased.
    field_list = [fdict['field'] for fdict in field_dict_list]

    # Get the verbose field names, which will be used as column headers.
    def _get_verbose(fdict):
        if 'verbose' in fdict:
            return fdict['verbose']
        else:
            return string.capwords(fdict['field'],'_').replace('_',' ')
    field_verbose_names = [_get_verbose(fdict) for fdict in field_dict_list]

    # A list of dicts containing the order of each column and the field titles
    # for each column, used for configuring jquery.datatables.js
    obj_field_config = [{
        'mData': name,
        'sTitle': verbose_name
    } for (name, verbose_name) in zip(field_list, field_verbose_names)]

    # Package the result.
    return json.dumps({
        OBJ_LIST: fe_obj_list,
        'field_config': obj_field_config
    })



def adapt_model_instance_to_frontend(model_instance, field_info={}, **kwargs):
    """Adapts a single model instance to the frontend representation.

    Args:
        model_instance: An instance of a Model object. The model class must
            implement a get_field_order() method.
        field_info: If called recursively from
            get_model_field_fe_representation(), function also be passed
            with field_info keys from parent_model.get_field_order(). This
            can decorate the serialized model with information like CSS class,
            state, instructions on how to render in datatable_component.js,
            etc.

    Returns:
        A dictionary representation of the model. May contained nested
            objects.
    """
    # The model class.
    model_type = type(model_instance)

    # The visible fields of the model.
    visible_field_names = [f['field']
            for f in model_type.get_field_order(**kwargs)]
    visible_field_dict = {f['field']: f
            for f in model_type.get_field_order(**kwargs)}

    # Get (key, value) pairs for visible fields.
    visible_field_pairs = [
            (field, get_model_field_fe_representation(
                    model_instance, field, visible_field_dict[field],
                    **kwargs))
            for field in visible_field_names]

    # Other values.
    other_pairs = []
    for key in HIDDEN_PAIR_FIELDS:
        if hasattr(model_instance, 'custom_getattr'):
            hidden_value = model_instance.custom_getattr(key)
            if hidden_value == UNDEFINED_STRING:
                continue
            other_pairs.append((key, hidden_value))
        elif hasattr(model_instance, key):
            other_pairs.append((key, getattr(model_instance, key)))

    # Add in keys from field_info, which are inherited from parent model, if
    # this function is called recursively from
    # get_model_field_fe_representation().
    if field_info:
        other_pairs.extend(field_info.items())

    # Wrap the results in a dictionary.
    return dict(visible_field_pairs + other_pairs)


def get_model_field_fe_representation(model_obj, field, field_info={},
        **kwargs):
    """Returns the best frontend representation for a model field that is
    implemented.

    This method allows recursively diving into models.
    """
    if hasattr(model_obj, 'custom_getattr'):
        model_field = model_obj.custom_getattr(field)
    else:
        model_field = getattr(model_obj, field)

    # Maybe special handling if ModelField is of special type.
    if isinstance(model_field, Model):
        return adapt_model_instance_to_frontend(model_field, field_info)
    elif model_field.__class__.__name__ == 'ManyRelatedManager':
        return [adapt_model_instance_to_frontend(m, field_info)
                for m in model_field.all()]
    elif isinstance(model_field, QuerySet):
        return [adapt_model_instance_to_frontend(m, field_info)
                for m in model_field]
    elif isinstance(model_field, datetime.datetime):
        return model_field.strftime("%Y-%m-%d %H:%M:%S")

    # Default. No further special handling needed.
    return str(model_field)
