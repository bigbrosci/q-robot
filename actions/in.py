#!/usr/bin/env python3
from incar import *
import sys 

dict_tasks, dict_task_groups = analyze_tasks(tasks)
generate_incar(standard_incar, dict_tasks, dict_task_groups)

