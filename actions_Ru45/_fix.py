with open('9_GA_MKM_triangle_sites.py', 'r') as f:
    lines = f.readlines()

# Find the line with "from cluster import *"
new_lines = []
inserted = False
for i, line in enumerate(lines):
    if 'from cluster import *' in line and not inserted:
        # Insert sys.path setup before this line
        new_lines.append("import sys\n")
        new_lines.append("sys.path.insert(0, '/home/qli/bin/q-robot/brain')\n")
        new_lines.append("\n")
        inserted = True
        new_lines.append(line)
    elif not (i < 8 and ('sys' in line or 'try' in line or 'except' in line or 'q_robot_path' in line or 'if q_robot' in line)):
        new_lines.append(line)
    elif i >= 8:
        new_lines.append(line)

with open('9_GA_MKM_triangle_sites.py', 'w') as f:
    f.writelines(new_lines)

print("Fixed!")
