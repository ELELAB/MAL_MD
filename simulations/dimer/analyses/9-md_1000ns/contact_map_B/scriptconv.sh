#!/usr/bin/env python

import sys

def unquote(s):
    return s[1+s.find('"'):s.rfind('"')]

def uncomment(s):
    return s[2+s.find('/*'):s.rfind('*/')]

def col(c):
    color = c.split('/*')
    value = unquote(color[1])
    color = unquote(color[0]).split()
    sys.stderr.write("%s: %s %s %s\n"%(c.strip(),color[0], color[1], value))
    return color[0], value

# Open the xpm file for reading
xpm = open(sys.argv[1])

# Read in lines until we fidn the start of the array
meta = [xpm.readline()]
while not meta[-1].startswith("static char *gromacs_xpm[]"):
    meta.append(xpm.readline())

# The next line will contain the dimensions of the array
dim = xpm.readline()
# There are four integers surrounded by quotes
nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]

# The next dim[2] lines contain the color definitions
# Each pixel is encoded by dim[3] bytes, and a comment
# at the end of the line contains the corresponding value
colors = dict([col(xpm.readline()) for i in range(nc)])

for i in xpm:
    if i.startswith("/*"):
        continue
    j = unquote(i)
    z = [colors[j[k:k+nb]] for k in range(0,nx,nb)]
    sys.stdout.write(" ".join(z)+"\n")

###
