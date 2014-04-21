#!/bin/bash

if [ -x "$PWD/manage.py" ]; then
  echo ok
else
  echo "manage.py not found. Please run this script from directory where "`
        `"manage.py is located."
  exit 1
fi

./manage.py test --settings=tests.integration_test_settings tests/integration
