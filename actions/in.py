#!/usr/bin/env python3
from incar import *
import sys 
from difflib import SequenceMatcher


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

tasks = [i.lower() for i in sys.argv[1:]]

tasks_not_recorded = []
for i in tasks:
    if i not in tasks_recorded:
        tasks_not_recorded.append(i)

if len(tasks_not_recorded) >= 1:
    for i in tasks_not_recorded:
        for j in tasks_recorded:
            if similar(i, j) >= 0.6:
                print('\n%s is not supported, perhaps you want to run task: %s' %(i, j))
    print('\nSupported tasks: %s\n' %(' '.join(tasks_recorded)))     
    print('Please use the right tasks above and run this command again.')
    exit()
        
dict_tasks, dict_task_groups = analyze_tasks(tasks)
generate_incar(standard_incar, dict_tasks, dict_task_groups)

