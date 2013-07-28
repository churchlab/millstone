"""Adapter - Converts django models to Front-end format.

    Calls on the model's get_field_order method to know which fields that 
    are available for display and what order. Returns a config json
    that dumps the field config and the json data for model/table.
"""

import json

from django.db.models import Model

def get_adapter(model, filters={}):
    """Converts django models to Front-end format.

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

    # Fetch all objects in this model and render it into the dom as json.
    # The data will be displayed to the user via the javascript DataTables
    # component.

    # Get all objects that pass the filter.
    obj_list = model.objects.filter(**filters)

    # Get a list of fields required for displaying the objects, in the order
    # in which they should be displayed.
    field_list = model.get_field_order()

    # Get the verbose field names, which will be used as column headers.
    field_verbose_names = [model._meta.get_field(field).verbose_name
            for field in field_list]

    def _get_fe_representation(model_obj, field):
        """Returns the best frontend representation for a model field that is
        implemented.

        For unicode fields, we can just return them as strings. However, it's
        more complicated for fields that are models themselves. In that case,
        we look for a get_fe_representation method on the model and use that if
        possible, else resort to casting the model to a string.
        """
        model_field = getattr(model_obj,field)
        if isinstance(model_field, Model) and hasattr(
                model_field, 'get_fe_representation'):
            return model_field.get_fe_representation()
        return str(model_field)

    # A list of dicts with object data, where each dict is one object
    # and all the fields required for front-end display.
    fe_obj_list = [dict([(field, _get_fe_representation(obj, field))
        for field in field_list] + [('href', 'http://google.com')])
            for obj in obj_list]

    # A list of dicts containing the order of each column and the field titles
    # for each column, used for configuring jquery.datatables.js
    obj_field_config = [{'mData': name, 'sTitle': verbose_name
        } for (name, verbose_name) in zip(field_list, field_verbose_names)]

    # The json result.
    obj_list_json = json.dumps({
        'obj_list': fe_obj_list,
        'field_config': obj_field_config
    })

    return obj_list_json
