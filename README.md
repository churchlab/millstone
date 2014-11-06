Millstone
==================

Millstone is a distributed bioinformatics software platform designed to facilitate genome engineering for synthetic biology. Automate iterative design, analysis, and debugging for projects involving hundreds of microbial genomes.

## Dependencies

* Python 2.7.3
* Perl 5 (for JBrowse)
* Java 1.7
* Postgresql 9.3 (only this version has been tested; on Mac we recommend [Postgres.app](http://postgresapp.com/))
* Unafold (http://dinamelt.rit.albany.edu/download.php)
* Python deps: See requirements.txt / instructions below
* RabbitMQ (not required to pass tests, see below)


## Installation

These are directions for installing Millstone on your own machine, and is **meant for advanced users who want a custom installation**. If you want to deploy a pre-configured Millstone instance to the cloud, then read [Getting Started with Millstone](https://github.com/churchlab/millstone/wiki/Getting-Started-with-Millstone) on the wiki.

**Before continuing, make sure all above dependencies are installed.** On Mac, we prefer [Homebrew](http://brew.sh/) for package management, and use it in the instructions below.

### Cloning the repository

Before installing, you must install `git` and clone the latest version of Millstone from GitHub. GitHub has [information on setting up git](https://help.github.com/articles/set-up-git#setting-up-git). Once Git is installed, you can clone the repository with:

        $ git clone https://github.com/churchlab/millstone.git <millstone_installation dir>
        $ cd <millstone_installation dir>

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
>If you want to use a version of python different from the OS default you can specify the python binary with the '-p' option:
>
>        $ virtualenv -p /usr/local/bin/python2.7 ~/pyenvs/genome-designer-env

3. Activate the environment in the shell. This will use `python` and other binaries like `pip` that are located your pyenv. You should do this whenever running any python/django scripts.

        $ source ~/pyenvs/genome-designer-env/bin/activate .

4. Install the dependencies in your virtual environment. We use the convention of running `pip freeze` to a .txt file containing a list of requirements.
Most users will want to do:

        $ pip install -r requirements/deploy.txt

If you plan on editing the code, you should run:

        $ pip install -r requirements/dev.txt

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

        $ cd jbrowse
        $ ./setup.sh
        $ cd ..


2. Install nginx if it's not already installed and copy
   or symlink the config file to nginx sites-enabled dir.

   > NOTE: On Mac, the sites-enabled dir is not present on the nginx version installed 
   > with brew, so the directory has to be added and included into the

   Unix:
        
        $ sed -i.orig "s:/path/to/millstone:$(pwd):g" config/jbrowse.local.nginx
        $ ln -s config/jbrowse.local.nginx /etc/nginx/sites-enabled

   Mac (run these commands from the project root):

        $ sed -i.orig "s:/path/to/millstone:$(pwd):g" config/jbrowse.local.nginx
        $ perl -pi.orig -e '$_ .= "    include /usr/local/etc/nginx/sites-enabled/*;\n" if /^http/' /usr/local/etc/nginx/nginx.conf
        $ sudo mkdir -p /usr/local/etc/nginx/sites-enabled
        $ sudo ln -s `pwd`/config/jbrowse.local.nginx /usr/local/etc/nginx/sites-enabled/millstone
        

3. Restart nginx.

   Unix:

        $ sudo service nginx restart

   Mac:

        $ ln -sfv /usr/local/opt/nginx/*.plist ~/Library/LaunchAgents
        $ launchctl load ~/Library/LaunchAgents/homebrew.mxcl.nginx.plist
        $ # to reload: sudo nginx -s reload
        
3. Install Perl local::lib module from CPAN for Jbrowse:

        $ sudo cpan install local::lib

4. Check that JBrowse is working locally by visiting:

      <http://localhost/jbrowse/index.html?data=sample_data/json/volvox>

*NOTE:* If upon running the Millstone application or its tests you observe errors
related to missing perl modules, you should also install them with `cpan`.

### Async Queue - RabbitMQ backend for Celery (optional for dev)

*NOTE:* Tests should pass without RabbitMQ setup so okay to skip this at first.

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

        $ brew install rabbitmq

    After install, you can run the server with:

        $ sudo /usr/local/sbin/rabbitmq-server

    Further Mac instructions are [here](http://www.rabbitmq.com/install-homebrew.html).
    
### Run the Millstone setup script.

The following installs various third-party bioinformatics tools and sets up
JBrowse.

    $ cd genome_designer
    $ ./millstone_setup.py

### Configuring PostgreSQL database for Millstone.

**NOTE:** If you make local changes, be sure to put them in a file called `genome_designer/conf/local_settings.py`.  You should not modify `global_settings.py`.

1. _(Mac Only)_ If you are using a fresh Postgres install, you may need to initialize the database:

        $ initdb /usr/local/var/postgres -E utf8
        $ pg_ctl -D /usr/local/var/postgres -l logfile start 

2. _(Mac Only)_ Since most new postgres Mac installations (both via brew and Postgres.app) do not have a `postgres` admin user, you will need to [modify your DATABASES variable](https://github.com/churchlab/millstone/issues/433) in `genome_designer/conf/local_settings.py`.

3. Navigate to the the `genome_designer/` dir.

4. Bootstrapping the database will automatically add the user, db, and permissions:

        
        $ python scripts/bootstrap_data.py -q

## Tests

We have two kinds of tests: unit and integration. Unit tests are intended
to be more specific tests of smaller pieces of the code, while integration
tests attempt to connect multiple pieces. Also, the integration tests actually
start celery worker intstances to simulate what happens in an async
environment while our unit tests use `CELERY_ALWAYS_EAGER = True` to mock out
testing celery.

We currently use [django-nose](https://pypi.python.org/pypi/django-nose) for
running tests, which provides a better interface than Django's native testing
setup (although this might not be true with the latest Django).

To run unit tests:

    (venv)$ ./scripts/run_unit_tests.sh

To run integration tests:

    (venv)$ ./scripts/run_integration_tests.sh

Nose also allows us to run tests only in specific modules.

In order to run only the tests in, say, the `main` app directory, run:

    (venv)$ ./scripts/run_unit_tests.sh main

For integration tests, we haven't figured out the optimal syntax in the test
script so to run individual tests, you'll need to do it this more explicit way:

    (venv)$ ./manage.py test --settings=tests.integration_test_settings tests/integration/test_pipeline_integration.py:TestAlignmentPipeline.test_run_pipeline

The same form works for unit tests, just use `--settings=tests.test_settings`

Note, in the following examples we use a standard `manage.py test` root, but you should adhere to the
examples above.

To run a single test module, run:

    (venv)$ python manage.py test main.tests.test_models

To run a single test case, e.g.:

    (venv)$ python manage.py test scripts/tests/test_alignment_pipeline.py:TestAlignmentPipeline.test_create_alignment_groups_and_start_alignments

To reuse the Postgresql database, wiping it rather than destroying and creating each time, use:

    (venv)$ REUSE_DB=1 ./manage.py test

Note that for some reason integration tests currently fail if run with the form:

    (venv)$ REUSE_DB=0 ./scripts/run_integration_tests.sh
    
Make sure you have R and unafold installed to avoid errors.

### Integration Tests

We recently introduced the concept of integration tests to our code. Previously,
many of our unit tests outgrew their unit-ness, but we were still treating them
like so.

We created an IntegrationTestSuiteRunner where the main difference is that
we start up a celery server that handles processing tasks. We are migrating
tests that should really be integration tests to be covered under this label.

When adding a test (see below), if your test touches multiple code units, it's
likely that's more appropriate to put it under integration test coverage. We'll
add notes shortly about how to add new integration tests.

To run integration tests, use this command. This uses nose so you can use
the same options and features as before.

    (venv)$ ./scripts/run_integration_tests.sh

HINT: When debugging integration tests, it may be necessary to manually clean
up previously stared `celerytestworker`s. There is a script to do this for you:

    $ ./scripts/kill_celerytestworkers.sh
 

### Debugging Tests / Dealing with Craziness

Our test framework isn't perfect. Here are some potential problems and other
hints that might help.

#### Celery workers not starting?

You might see this error:

    AssertionError: No running Celery workers were found.

So far, we're aware of a couple reasons you might see this:

* Another integration test in the same run failed, causing all
  subsequent integration tests to fail.
* celerytestworkers were not shut down (use script above to kill them).

#### `ps aux` and `grep` are your friends

To see running celery processes:

    ps aux | grep celery

To see running integration test:

    ps aux | grep python.*integration

To kill the process associated with the integration test (e.g. 777)

    kill 777

To kill orphaned celerytestworker processes, we actually have a script:

    ./scripts/kill_celerytestworkers.sh

#### Running individual tests

This command runs a specific integration test and doesn't capture stdout:

    ./manage.py test -s --settings=tests.integration_test_settings tests/integration/test_pipeline_integration.py:TestAlignmentPipeline.test_run_pipeline


### Adding Tests

(Right now, this documentation is only for unit tests. Information for
integration tests is coming soon.)

Nose automatically discovers files with names of the form `test_*.py` as test
files.

## Running the application

0. Activate your virtualenv, e.g.:

        $ source ~/pyenvs/genome-designer-env/bin/activate .

1. Navigate to the the `genome_designer/` dir.

2. From one terminal, start the celery server.

        (venv)$ ./scripts/run_celery.sh

3. Open another terminal and start the django server.

        (venv)$ python manage.py runserver

4. Visit the url <http://localhost:8000/> to see the demo.

## Bootstrapping Test Data

First make sure Celery is running. In another terminal do:

    (venv)$ ./scripts/run_celery.sh

From the `genome_designer` directory, run:

    (venv)$ python scripts/bootstrap_data.py full

NOTE: This will delete the entire dev database and re-create it with the
hard-coded test models only.

## Debugging

### Accessing the Postgresql database

On Ubuntu, if your database is called `gdv2db`:

    sudo -u postgres psql gdv2db
    
To debug tests with pdb, add `pdb.set_trace()` checkpoints and use a command similar to:

    REUSE_DB=1 ./manage.py test -s --pdb --pdb-failures main/tests/test_xhr_handlers.py:TestGetVariantList.test__basic_function

### Profiling code

The `debug.profiler` module contains a `profile` decorator that can be added to a function. For example, to debug a view:

1. Update local_settings.py with your profiler logs destination by setting the PROFILE_LOG_BASE attribute, e.g.:

        PROFILE_LOG_BASE = '/path/to/logs'

2. Add @profile('log_file_name') in front of the method you want to profile, e.g.:

        @profile('mylog')
        def my_view(request):
            ...

3. Use the `debug/inspect_profiler_data.py` convenience script to parse the data, e.g.:

        python inspect_profiler_data.py /path/to/log/mylog
