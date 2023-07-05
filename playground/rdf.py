#! /usr/bin/env python3
import sys
import json
import pprint
import rdflib
# ATTN: on Ubuntu, Centos etc should be run python3 semweb.py !!

myfile = sys.argv[1]
g = rdflib.Graph()
# default format: 'xml'
# from https://media.readthedocs.org/pdf/rdflib3/latest/rdflib3.pdf
# the formats recognized by the serializers:
# n3, nquads, nt, turtle
# trig, trix
# xml, pretty-xml

# formats successfilly tested: nt, n3, turtle
# no luck with nquads

# parsers accept extra formats compared to serializers, e.g. rdfa, microdata, etc
result = g.parse(myfile, format='n3')
print("graph has %s statements." % len(g))

for subj, pred, obj in g:
   if(subj, pred, obj) not in g:
       raise Exception("It better be!")

s = g.serialize(destination=myfile, format='turtle')
# desitnation is optional
# if provided the return value is 'None', the stream sent streight into the file
print(s)
