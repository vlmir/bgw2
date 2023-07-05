#! /usr/bin/env python3
import sys
import json
import pprint
import rdflib
import networkx as nx
import pronto
import owlready2
import obonet

# ATTN: on Ubuntu, Centos etc should be run python3 semweb.py !!

###############################################################################
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
myfile = '/home/mironov/git/usr/bgw/tdata/obof/tgo.ttl'
result = g.parse(myfile, format='n3')
print("graph has %s statements." % len(g))

for subj, pred, obj in g:
   if(subj, pred, obj) not in g:
       raise Exception("It better be!")

s = g.serialize(destination='webtest.ttl', format='turtle')
# desitnation is optional
# if provided the return value is 'None', the stream sent streight into the file

###############################################################################

myfile = '/home/mironov/git/usr/bgw/data/onto/ro.obo'
myfile = '/home/mironov/git/usr/bgw/data/onto/go-basic.obo'
myfile = '/home/mironov/git/usr/bgw/tdata/obof/tgo.obo'
onto = pronto.Ontology(myfile)
#onto = pronto.Ontology('http://purl.obolibrary.org/obo/mondo.obo')
# does not import Typedefs !!
print(len(onto))
with open('webtest.obo', 'w') as fp:
    fp.write(onto.obo)
#for r in Relationship.topdown():
#    print(r)
if 'GO:0008150' in onto:
    print('True')
else:
    print('False')

print(onto['GO:0008150'].relations)

###############################################################################

graph = obonet.read_obo(myfile)
print(len(graph))
print(graph.number_of_edges())
if nx.is_directed_acyclic_graph(graph):
    print('True')
else:
    print('False')
trmid = 'GO:0008150' 
if trmid in onto:
    print('True')
else:
    print('False')
pprint.pprint(graph)
print(nx.descendants(graph, trmid))
print(nx.ancestors(graph, trmid))
data = nx.node_link_data(graph)
data = nx.adjacency_data(graph)
#data = nx.jit_data(graph) # not implemented for multigraph type
#data = nx.tree_data(graph, trmid) # TypeError: G is not a tree.
s = json.dumps(data)
with open('test.json', 'w') as fp:
    fp.write(s)
