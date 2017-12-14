#!/usr/bin/env python2.7
# r2 .. just takes a bed file wkexnrt.py a file to WorK on exonerate output
# operates on output of the exonerate program.
# Rather simply I'm afraid, just regex'ing the Target range line
import sys, regex
from collections import defaultdict

def main():
    """ Generate a table from mash results """
    argquan=len(sys.argv)
    if argquan != 2:
        print "This script requires one argument: a file listing of the mash re sults files"
        sys.exit(2)

    # regex
    # Target range: 138 -> 222
    RGX2=regex.compile(r'^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+).*')
    with open(sys.argv[1]) as x: fl = x.read().splitlines()
    D=defaultdict(list) # this method from collections module, allows undefined keys to be initialiased cleanly.
    L=[]
    for i in fl:
        m2=RGX2.match(i)
        mm0= m2.groups(1)[4]
        if(m2):
            pd0= int(m2.groups(1)[1]) # come out as strings of course ... pd, pair distance ... this we can convert to integer.
            pd2= int(m2.groups(1)[2]) 
            if pd2 > pd0:
                d= pd2 - pd0
            else:
                d = -1 * (pd2 -pd0)
            D[mm0]=(pd0, pd2, d) # so each sequence ID will be a dict key value with a list of 2-uples.


    DSZ=len(D)
    for k,v in D.items():
        print "%s len = %i" % (k, D[k][2])

if __name__=='__main__':
    main()
