from django.http import Http404
from django.http import HttpResponse
from django.http import HttpResponseBadRequest
from django.http import HttpResponseRedirect
from django.shortcuts import get_object_or_404
from django.shortcuts import render
from django.views.decorators.http import require_POST
from django.conf import settings
import psycopg2
import json
import os
import kombu


def get_local_settings_path():
    path = os.path.join(settings.BASE_DIR, "../../genome_designer/conf/local_settings.py")
    assert os.path.isfile(path)
    return path

def home(request):
    path = get_local_settings_path()
    with open(path) as f:
        code = compile(f.read(), path, 'exec')
    config = {}
    exec code in config

    c = {
        'dbname': config.get('DATABASES', {}).get('default', {}).get('NAME', 'genome-designer'),
        'user': config.get('DATABASES', {}).get('default', {}).get('USER', 'genome-designer'),
        'password': config.get('DATABASES', {}).get('default', {}).get('PASSWORD', ''),
        'host': config.get('DATABASES', {}).get('default', {}).get('HOST', ''),
    }
    return render(request, 'base.html', c)

@require_POST
def test(request):
    user = request.POST['user']
    password = request.POST['password']
    host = request.POST['host']
    dbname = request.POST['dbname']

    try:
        connstring = "dbname='%s' user='%s' host='%s' password='%s' port=5432 connect_timeout=5" % (
            dbname, user, host, password)
        conn = psycopg2.connect(connstring)
    except Exception as e:
        print e
        response = {
            'status': 'error',
            'message': "Cannot connect to PostgresSQL: " + str(e),
            'connstring': connstring
        }
        return HttpResponse(json.dumps(response),
            content_type='application/json')

    try:
        connstring = 'amqp://%s:%s@%s:5672//' % (user, password, host)
        conn = kombu.Connection(connstring)
        conn.connect()
    except Exception as e:
        print e
        response = {
            'status': 'error',
            'message': "Cannot connect to RabbitMQ: " + str(e),
            'connstring': connstring
        }
        return HttpResponse(json.dumps(response),
            content_type='application/json')

    response = {
        'status': 'success',
    }
    return HttpResponse(json.dumps(response),
        content_type='application/json')

def save(request):
    user = request.POST['user']
    password = request.POST['password']
    host = request.POST['host']
    dbname = request.POST['dbname']

    db = {
        'default': {
            'ENGINE': 'django.db.backends.postgresql_psycopg2',
            'NAME': dbname,
            'USER': user,
            'PASSWORD': password,
            'HOST': host,
            'PORT': 5432,
        }
    }
    rabbitmq = "amqp://%s:%s@%s:5672//" % (user, password, host)
    config = "DATABASES = %s \n" % repr(db)
    config += "BROKER_URL = %s \n" % repr(rabbitmq)
    path = get_local_settings_path()

    with open(path, "w") as f:
        f.write(config)

    os.system("supervisorctl restart celery")

    response = {
        'status': 'success'
    }
    return HttpResponse(json.dumps(response),
        content_type='application/json')
