genome-designer-v2
==================

Second generation of Genome Designer

## Dependencies

* Postgresql 9.3 (only this version has been tested)
* R (for Picard)
* Python deps: See requirements.txt / instructions below

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

NOTE: Watch changes to requirements.txt and re-run the install command when
collaborators add new dependencies.

### Git Submodules

We currently submodule `jbrowse` and perhaps will do so with other tools in the
future. Specifically, we have submoduled a forked copy of `jbrowse` at a
specific commit. To checkout the appropriate submodule states, run:

    $ git submodule update --init --recursive

This will pull the submodules and also pull any of their submodules.


### JBrowse (continued)

After installing JBrowse via the Git submodule route described above, you
need to do the following to get JBrowse up and running:

*NOTE:* Only step 1 is necessary to get the tests to pass. The later steps need to be updated.

1. Run the JBrowse setup script inside the `jbrowse` dir.

        $ ./setup.sh


2. Install nginx if it's not already installed (use brew on Mac OSX) and copy
   or symlink the config file to nginx sites-enabled dir.

   Unix:

        $ ln -s config/jbrowse.local.nginx /etc/nginx/sites-enabled

   Mac:

        $ sudo mkdir -p /usr/local/etc/nginx/sites-enabled
        $ sudo ln -s config/jbrowse.local.nginx /usr/local/etc/nginx/sites-enabled

3. Restart nginx.

   Unix:

        $ sudo service nginx restart

   Mac:

        $ ln -sfv /usr/local/opt/nginx/*.plist ~/Library/LaunchAgents
        $ launchctl load ~/Library/LaunchAgents/homebrew.mxcl.nginx.plist

4. Check that JBrowse is working locally by visiting:

      <http://localhost/jbrowse/index.html?data=sample_data/json/volvox>

If upon running the Genome Designer application or its tests you observe errors
related to missing perl modules, you can install them with `cpanm`, e.g.:

        $ cpanm local::lib


### Async Queue - RabbitMQ backend for Celery (optional for dev)

*NOTE:* Tests sould pass without RabbitMQ setup so okay to skip this at first.

Asynchronous processing is necessary for many of the analysis tasks in this
application.  We use the open source project celery since it is being actively
developed and has a library for integrating with Django. Celery requires a
message broker, for which we use RabbitMQ which is the default for Celery.

1. Install Celery

    The `celery` and `django-celery` packages are listed in
    requirements.txt and should be installed in your virtualenv following the
    instructions above.

2. Install RabbitMQ - On Ubuntu, install using sudo:

        $ sudo apt-get install rabbitmq-server

    Full instructions are [here](http://www.rabbitmq.com/download.html).

    On Mac, homebrew can be used:

        $ sudo brew install rabbitmq

    After install, you can run the server with:

        $ sudo /usr/local/sbin/rabbitmq-server

    Further Mac instructions are [here](http://www.rabbitmq.com/install-homebrew.html).



### Other third-party tools

For most third-party tools that the application depends on, we've put them in a
Dropbox, as either Mac OS X or Linux binaries. Run `setup.py` which should
create a `tools/` dir and will download the correct build for your system.
We recommend using our builds, as they have been tested with the rest of
the application.


## Running the application

0. Activate your virtualenv, e.g.:

        $ source ~/pyenvs/genome-designer-env/bin/activate .

1. Navigate to the the `genome_designer/` dir.

2. From one terminal, start the celery server.

        (venv)$ ./run_celery.sh

3. Open another terminal and start the django server.

        (venv)$ python manage.py runserver

4. Visit the url <http://localhost:8000/> to see the demo.


## Configuring PostgreSQL database.

See instructions for setting up PostgresSQL database:
<https://overlappingminds.com/thoughts/069accf7-ccb4-4b7a-bd40-022c49a053cd>

*NOTE:* Be sure to make local db config changes in `conf/local_settings.py`.

In order to run tests with Postgres, your user will need CREATE permissions.
Otherwise you might get an error creating a database.
You can grant these by logging into the Posgres shell and running:

    ALTER USER django CREATEDB;


## Tests

We currently use [django-nose](https://pypi.python.org/pypi/django-nose) for
testing. This package should be seamlessly hooked up to Django's normal testing
so you can do the standard `manage.py` command:

    (venv)$ python manage.py test

Nose also allows us to run tests only in specific modules.

In order to run only the tests in, say, the `main` app directory, run:

    (venv)$ python manage.py test main

And for only the tests in `scripts` call:

    (venv)$ python manage.py test scripts

To run a single test module, run:

    (venv)$ python manage.py test main.tests.test_models

To run a single test case, e.g.:

    (venv)$ python manage.py test scripts/tests/test_alignment_pipeline.py:TestAlignmentPipeline.test_create_alignment_groups_and_start_alignments

To reuse the Postgresql database, wiping it rather than destroying and creating each time, use:

    (venv)$ REUSE_DB=1 ./manage.py test


### Adding Tests

Nose automatically discovers files with names of the form `test_*.py` as test
files.

## Bootstrapping Test Data

From the `genome_designer` directory, run:

    (venv)$ python scripts/bootstrap_data.py

NOTE: This will delete the entire dev database and re-create it with the
hard-coded test models only.


## Profiling code

The `debug.profiler` module contains a `profile` decorator that can be added to a function. For example, to debug a view:

1. Update local_settings.py with your profiler logs destination by setting the PROFILE_LOG_BASE attribute, e.g.:

        PROFILE_LOG_BASE = '/path/to/logs'

2. Add @profile('log_file_name') in front of the method you want to profile, e.g.:

        @profile('mylog')
        def my_view(request):
            ...

3. Use the `debug/inspect_profiler_data.py` convenience script to parse the data, e.g.:

        python inspect_profiler_data.py /path/to/log/mylog
