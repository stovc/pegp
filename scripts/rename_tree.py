"""Fix the naming issue introduced by MEGA.
Removes the introduced `/...` part for each tip
"""

from ete3 import Tree

t = Tree("tree-nj.nwk")

for l in t.iter_leaves():
    l.name = l.name.split('/')[0]

out = open('tree-nj-rn.nwk', 'w')

out.write((t.write()))
out.close()
