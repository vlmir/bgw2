#! /usr/bin/env python3
import sys
from subprocess import call

path = sys.argv[1]
stream = open(path)
with open(path) as fp:
#    url = fp.readline()
#    while url:
#        call(['wget', url.strip()])
#        url = fp.readline()
## alternatively and BETTER:
#    for count, url in enumerate(fp):
#        call(['wget', url.strip()])

## yet another way, from: https://docs.python.org/3.3/tutorial/inputoutput.html
    for url in fp:
        call(['wget', url.strip()]) # the best!
