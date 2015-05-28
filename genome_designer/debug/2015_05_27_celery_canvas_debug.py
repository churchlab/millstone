"""Celery failure case.
"""

from celery import Celery
from celery import group

app = Celery('2015_05_27_celery_canvas_debug')
app.config_from_object('debug_celeryconfig')


@app.task
def t(v):
    print v
    return v


if __name__ == '__main__':

    ### From Github issue

    # g = group(
    #     t.si(1),
    #     chain(
    #         t.si(2),
    #         t.si(3)
    #     )
    # )

    # print 'Now get results'
    # print chain(g).delay().get()  # this works
    # # chain(g, t.si(4)).delay().get()  # this works
    # # chain(t.si(0), g).delay().get()  # this fails

    ### Reproducing Millstone bug.

    # This works.
    group_1 = group([
        t.si('group_1_task_1'),
        t.si('group_1_task_2'),
    ])
    group_2 = group([
        t.si('group_2_task_1'),
    ])
    addtl_task_1 = t.si('addtl_task_1')
    addtl_task_2 = t.si('addtl_task_2')
    pipeline = group_1 | group_2 | addtl_task_1 | addtl_task_2
    assert pipeline.delay().get() == 'addtl_task_2'

    # This fails.
    group_1 = group([
        t.si('group_1_task_1'),
    ])
    group_2 = group([
        t.si('group_2_task_1'),
        t.si('group_2_task_2'),
    ])
    addtl_task_1 = t.si('addtl_task_1')
    addtl_task_2 = t.si('addtl_task_2')
    pipeline = group_1 | group_2 | addtl_task_1 | addtl_task_2
    assert pipeline.delay().get() == 'addtl_task_2'

    # This also fails.
    group_1 = group([
        t.si('group_1_task_1'),
        t.si('group_1_task_2')
    ])
    group_2 = group([
        t.si('group_2_task_1'),
        t.si('group_2_task_2') # NEW
    ])
    addtl_task_1 = t.si('addtl_task_1')
    addtl_task_2 = t.si('addtl_task_2')
    pipeline = group_1 | group_2 | addtl_task_1 | addtl_task_2
    assert pipeline.delay().get() == 'addtl_task_2'

    # Yet another fail case.
    # Apparently having a group as a header of any chord is broken.
    group_2 = group([
        t.si('group_2_task_1'),
        t.si('group_2_task_2'),
    ])
    addtl_task_1 = t.si('addtl_task_1')
    addtl_task_2 = t.si('addtl_task_2')
    pipeline = group_2 | addtl_task_1 | addtl_task_2
    assert pipeline().get() == 'addtl_task_2', pipeline().get()
