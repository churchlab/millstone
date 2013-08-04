"""Adapter - Converts django models to Front-end format.

    Calls on the model's get_field_order method to know which fields that 
    are available for display and what order. Returns a config json
    that dumps the field config and the json data for model/table.
"""

import json
import string

from django.db.models import Model, ManyToManyField

#from django.db.models.fields.related import ManyRelatedManager

def adapt_model_to_frontend(model, filters={}):
    """Converts django models to frontend format.

    Calls on the model's get_field_order method to know which fields that
    are available for display and what order.

    Args:
        model: The django model that we are adapting
        **kwargs: The filter conditions on the model (i.e. select
            which instances of the model to select.)

    Returns:
        A config json that dumps the field config and the json data for
        model/table display.
    """

    # Fetch all objects in this model and package as json to be rendered into
    # the dom. The data will be displayed to the user, e.g. via the javascript
    # DataTables component.

    # Get all objects that pass the filter.
    obj_list = model.objects.filter(**filters)

    # A list of dicts with object data, where each dict is one object
    # and all the fields required for front-end display.
    fe_obj_list = [adapt_model_instance_to_frontend(obj) for obj in obj_list]

    # Get a list of fields required for displaying the objects, in the order
    # in which they should be displayed.
    field_dict_list = model.get_field_order()
    
    # Each field is a dict with two keys, 'field' for field name and 'verbose'
    # for display name. Get each. If 'verbose' is missing, then make verbose 
    # be the field with _'s turned to spaces and Title Cased. 
    field_list = [fdict['field'] for fdict in field_dict_list]

    # Get the verbose field names, which will be used as column headers.
    get_verbose= lambda fdict: (fdict['verbose'] if 'verbose' in fdict else 
        string.capwords(fdict['field'],'_').replace('_',' '))    
    field_verbose_names = [get_verbose(fdict) for fdict in field_dict_list]


    # A list of dicts containing the order of each column and the field titles
    # for each column, used for configuring jquery.datatables.js
    obj_field_config = [{
        'mData': name,
        'sTitle': verbose_name
    } for (name, verbose_name) in zip(field_list, field_verbose_names)]


    # Package the result.
    return json.dumps({
        'obj_list': fe_obj_list,
        'field_config': obj_field_config
    })


def adapt_model_instance_to_frontend(model_instance, field_info={}):
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
    visible_field_list = [f['field'] for f in model_type.get_field_order()]
    visible_field_dict = {f['field'] : f for f in model_type.get_field_order()}


    # Get (key, value) pairs for visible fields.
    visible_field_pairs = [(field, get_model_field_fe_representation(
        model_instance, field, visible_field_dict[field])) for field 
        in visible_field_list]

    # Other values.
    other_pairs = []
    if hasattr(model_instance, 'get_href'):
        other_pairs.append(('href', model_instance.get_href()))
    if hasattr(model_instance, 'uid'):
        other_pairs.append(('uid', model_instance.uid))
        
    # Add in keys from field_info, which are inherited from parent model, if
    # this function is called recursively from 
    # get_model_field_fe_representation().
    if field_info: 
        other_pairs.extend(field_info.items())
    
    print field_info
    print visible_field_pairs
    print other_pairs
    
    # Wrap the results in a dictionary.
    return dict(visible_field_pairs + other_pairs)


def get_model_field_fe_representation(model_obj, field, field_info={}):
    """Returns the best frontend representation for a model field that is
    implemented.

    This method allows recursively diving into models.
    
    """
    model_field = getattr(model_obj,field)
    
    if isinstance(model_field, Model):
        return adapt_model_instance_to_frontend(model_field, field_info)
    elif model_field.__class__.__name__ ==  'ManyRelatedManager':
        return [
            adapt_model_instance_to_frontend(m, field_info) for m in model_field.all()]
    return str(model_field)
