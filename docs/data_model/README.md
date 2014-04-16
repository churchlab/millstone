# Creating Data Models

We use the [graph_models](http://django-extensions.readthedocs.org/en/latest/graph_models.html) tool inside of [django-extensions](https://github.com/django-extensions/django-extensions).

## Installation

Make sure the library is installed:

    (venv) $ pip install django-extensions

Update `settings.py` to install the app:

    INSTALLED_APPS = (
        ...
        'django_extensions',
        ...
    )

## Workflow

To generate an initial graphviz dot file, run:

    (venv) $ ./manage.py graph_models main > millstonedb.dot

Convert to png (or whatever):

    $ dot -Tpng millstonedb.dot -o millstonedb.png

After observing, edit the .dot file by handle to customize.
