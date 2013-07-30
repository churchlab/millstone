"""Adapter - Converts django models to Front-end format.

    Calls on the model's get_field_order method to know which fields that 
    are available for display and what order. Returns a config json
    that dumps the field config and the json data for model/table.
"""

import json

from django.db.models import Model

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
    field_list = model.get_field_order()

    # Get the verbose field names, which will be used as column headers.
    field_verbose_names = [model._meta.get_field(field).verbose_name
            for field in field_list]

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


def adapt_model_instance_to_frontend(model_instance):
    """Adapts a single model instance to the frontend representation.

    Args:
        model_instance: An instance of a Model object. The model class must
            implement a get_field_order() method.

    Returns:
        A dictionary representation of the model. Maybe contained nested
            objects.
    """
    # The model class.
    model_type = type(model_instance)

    # The visible fields of the model.
    visible_field_list = model_type.get_field_order()

    # Get (key, value) pairs for visible fields.
    visible_field_pairs = [(field, get_model_field_fe_representation(
        model_instance, field)) for field in visible_field_list]

    # Other values.
    other_pairs = []
    if hasattr(model_instance, 'get_href'):
        other_pairs.append(('href', model_instance.get_href()))
    if hasattr(model_instance, 'uid'):
        other_pairs.append(('uid', model_instance.uid))

    # Wrap the results in a dictionary.
    return dict(visible_field_pairs + other_pairs)


def get_model_field_fe_representation(model_obj, field):
        """Returns the best frontend representation for a model field that is
        implemented.

        This method allows recursively diving into models.
        """
        model_field = getattr(model_obj,field)
        if isinstance(model_field, Model):
            return adapt_model_instance_to_frontend(model_field)
        return str(model_field)
