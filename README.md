genome-designer-v2
==================

Second generation of Genome Designer


## Installation

### Python environment

We recommend using [virtualenv](http://pypi.python.org/pypi/virtualenv) for
creating and managing a sandbox python environment. This strategy makes it easy
to stay up with requirements. Our requirements are listed in `requirements.txt`.
Follow the instructions below to setup your virtualenv and install the required
packages.

#### Setting up a virtual python environment

1. Install [virtualenv](http://www.virtualenv.org/en/latest/index.html) if you don't have it yet. (You may want to install [pip](http://pypi.python.org/pypi/pip/) first.)

2. Create a new virtual environment for this project. This virtual environment isn't part of the project so just put it somewhere on your machine.  I keep all of my virtual environments in the directory ~/pyenvs/.

        $ virtualenv ~/pyenvs/genome-designer-env

    If you want to use a version of python different from the OS default you can specify the python binary with the '-p' option:

        $ virtualenv -p /usr/local/bin/python2.7 ~/pyenvs/genome-designer-env

3. Activate the environment in the shell. This will use `python` and other binaries like `pip` that are located your pyenv. You should do this whenever running any python/django scripts.

        $ source ~/pyenvs/genome-designer-env/bin/activate .

4. Install the dependencies in your virtual environment. We've exported the requirements in the requirements.txt file. In theory, these should all be installable with the single command:

        $ pip install -r requirements.txt

However, in reality, this doesn't seem to work perfectly. In particular, it may
be necessary to install specific packages first.


## Running the Django server

0. Activate your virtualenv, e.g.:

        $ source ~/pyenvs/genome-designer-env/bin/activate .

1. Navigate to the the `genome_designer/` directory.

        $ python manage.py runserver


2. Visit the url <http://localhost:8000/> to see the demo.
