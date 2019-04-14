# Calculates the cumulative percent of principal components.
# So if you select 
# neig = 10
# you will get file with 10 numbers.
# They say how much does the 1st principal component contributes to the total movement.
# Then, how much do the 1st and the 2nd principal component together contribute to the total movement.
# Then, how much do the 1st, the 2nd and the 3rd principal component together contribute to the total movement...

# Usage: get_variance.py eigenval.xvg > values.dat

#!/usr/bin/python
import sys
import os
import re
#from Numeric import array, cumsum
import scipy
from numpy import *
neig = 10

def get_rmsf(rmsfname):
    fin = open(rmsfname)
    lines = fin.readlines()
    fin.close()
    datavec = []
    for line in lines:
        if line[0] == '#' or line[0] == '@':
            continue
        datavec.append(float(line.split()[1]))
    return datavec

def printfile(filename):
    earray = array(get_rmsf(filename))
    earrayperc = cumsum(earray)/sum(earray)*100
    print filename,
    j=1
    for i in earrayperc[0:neig]:
        print str(j)+"\t%8.3f" % i
        j=j+1
def printdir(dirname):
    patxvg = re.compile("xvg")
    filenames = os.listdir(dirname)
    for filename in filenames:
        if patxvg.search(filename): 
            printfile(dirname+'/'+filename)

if __name__ == '__main__':
    if not (len(sys.argv) == 2 or len(sys.argv) == 3):
        print "Usage: %s <xvgfilename> [dir]" % sys.argv[0]
        sys.exit()
    if len(sys.argv) == 2:
        print_perc_first_ten = printfile
    else:
        print_perc_first_ten = printdir
    print_perc_first_ten(sys.argv[1])
