#! /usr/bin/env python3
# TODO strip() after each split()
import sys
import json
import os.path
import csv
import pprint
DQ = '"'
END = ' .\n'
US = '--'
AS = '/' # attribute separator in URIs, not used yet
SS = '!' # slash substitute in rdf:Statement URIs
DS = ' + ' # separator for concatenating defs, gene.py and prot.py
KS = '+' # key separator, only up2mim.py and zeno.py
NS = '-' # not used yet

def xmap(pth2src, kind, vind, dlm='\t'):
    mapA = {
    }
    with open(pth2src) as fp:
        for line in fp:
            line = line.rstrip()
            if not line:
                continue
            skip = ['#', '!']
            null = ['-', '']
            if line[0] in skip:
                continue
            row = line.split(dlm)
            (key, val) = (row[kind], row[vind])
            if val in null:
                continue # TODO to be tested
            if(key in mapA.keys()):
                mapA[key][val] = 1 # sic! '+=' throws errors
            else:
                mapA[key] = {}
                mapA[key][val] = 1
    return(mapA)

def xfixlen(pth2src, fields, mydlm=':'):
# TODO re-test; generalize
    data = {}
    with open(pth2src) as fp:
        for line in fp:
            upac = ''
            mim = ''
            mod = ''
            for key in fields.keys():
                (ind1, ind2) = (fields[key][0], fields[key][1])
                val = line[ind1:ind2]
                if key == 'uniprot':
                    upac = val.replace(' ', '')
                if key == 'omim':
                    bits = val.split(mydlm)
                    mim = bits[-1][:-1]
                if key == 'change':
                    mod = val.replace(' ', '')
            mykey = upac + KS + mim
            if(mykey in data.keys()):
                data[mykey]['change'][mod] = 1
            else:
                data[mykey] = {}
                data[mykey]['change'] = {}
                data[mykey]['change'][mod] = 1
    return(data)

## TODO fix :-separated ids
#Em
#HGNC
#HostDB
#hsa
#MetaCyc
#SLP
#Smp_009580.1
#Smp_196150.1
#WUGSC
def idmap(pth2src, kind, vind, fields, mapA={}):
    mapkeys = sorted(mapA.keys())
    sind=1
    data = {
    }
    dlm = '\t'
    skip = ['#', '!']
    null = ['-', '']
    with open(pth2src) as fp:
        for line in fp:
            if line[0] in skip:
                continue
            line = line.rstrip()
            if not line: # an empty line
                continue
            row = line.split(dlm)
            (key, src, val) = (row[kind], row[sind], row[vind])
            if val in null:
                continue
            if key in null:
                continue
            if src in fields.keys():
                ns = fields[src]
            else:
                continue
            if mapkeys:
                if key not in mapkeys:
                    continue
            if key not in data.keys():
                data[key] = {}
            if ns not in data[key].keys():
                data[key][ns] = {}
            data[key][ns][val] = 1
    return(data)

def x2D(pth2src, dlm, mydlm, srckeys, myind, mapA={}):
    data = {
    }
    skip = ['#', '!']
    null = ['-', '']
    mykeys = sorted(mapA.keys())
    if not os.path.exists(pth2src): # TODO implement in all functions
        print('zeno.x2D: ERROR: not found file: ' + pth2src)
        return(0)
    with open(pth2src) as fp:
        for line in fp:
            if line[0] in skip:
                continue
            line = line.rstrip()
            if not line:
                continue
            row = line.split(dlm)
            indkeys = row[myind]

            for indkey in indkeys.split(mydlm):
                if indkey in null:
                    continue
                if len(mykeys):
                    if indkey not in mykeys:
                        continue
                if(indkey not in data.keys()):
                    data[indkey] = {}

                for ind in srckeys.keys():
                    field = srckeys[ind]
                    if (field not in data[indkey].keys()):
                        data[indkey][field] = {}
                    try:
                        vals = row[ind] # sic
                    except IndexError:
                        continue
                    if(vals): # sic
                        indkeys = vals.split(mydlm)
                        for val in indkeys:
                            if val in null:
                                continue
                            data[indkey][field][val] = 1
    return(data)

def xflex(pth2src, clmndlm, srckeys, myind, mapA={}):
    data = {
    }
    skip = ['#', '!']
    null = ['-', '']
    mapkeys = sorted(mapA.keys())
    if not os.path.exists(pth2src): # TODO implement in all functions
        print('zeno.xflex: ERROR: not found file: ' + pth2src)
        return(0)
    mydic = srckeys[myind]
    myns = sorted(mydic.keys())[0] # the only one
    mydlm = mydic[myns]
    with open(pth2src) as fp:
        for line in fp:
            if line[0] in skip:
                continue
            line = line.rstrip()
            if not line:
                continue
            chunks = line.split(clmndlm)
            mychunk = chunks[myind]

            # TODO accomodate for no delimiter
            for myval in mychunk.split(mydlm):
                myval = myval.strip()
                if myval in null:
                    continue
                if len(mapkeys):
                    if myval not in mapkeys:
                        continue
                if(myval not in data.keys()):
                    data[myval] = {}

                for ind in srckeys.keys():
                    if ind == myind:
                        continue
                    onedic = srckeys[ind]
                    onens = sorted(onedic.keys())[0]
                    onedlm = onedic[onens]
                    if (onens not in data[myval].keys()):
                        data[myval][onens] = {}
                    try:
                        chunk = chunks[ind] # sic
                    except IndexError:
                        continue
                    if not chunk: # sic
                        continue
                    if not onedlm: # e.g. ""
                        bit = chunk
                        if bit in null:
                            continue
                        data[myval][onens][bit] = 1
                    else:
                        bits = chunk.split(onedlm)
                        for bit in bits:
                            bit = bit.strip()
                            if bit in null:
                                continue
                            data[myval][onens][bit] = 1
    return(data)

def x2Dx2D(pth2src, dlm, mydlm, fields, ind1, ind2, mapA={}):
    data = {
    }
    mykeys = sorted(mapA.keys())
    with open(pth2src) as fp:
        for line in fp:
            skip = ['#', '!']
            null = ['-', '']
            if line[0] in skip:
                continue
            line = line.rstrip()
            if not line:
                continue
            row = line.split(dlm)

            key1s = row[ind1]
            key2s = row[ind2]
            for key1 in key1s.split(mydlm):
                if len(mykeys):
                    if key1 not in mykeys:
                        continue
                for key2 in key2s.split(mydlm):
                    if len(mykeys):
                        if key2 not in mykeys:
                            continue
                    if(key1 not in data.keys()):
                        data[key1] = {}
                    if(key2 not in data[key1].keys()):
                        data[key1][key2] = {}
                    for ind in fields.keys():
                        field = fields[ind]
                        if (field not in data[key1][key2].keys()):
                            data[key1][key2][field] = {}
                        vals = row[ind] # sic
                        if(vals): # sic
                            bits = vals.split(mydlm)
                            for val in bits:
                                if val in null:
                                    continue # TODO to be tested
                                data[key1][key2][field][val] = 1
    return(data)

def header(itemUs, oppys, appys, prns):
# itemUs returned via arg
    buff = ''
    for key in ('ObjectProperty', 'AnnotationProperty', 'Class'):
        if key == 'Class':
            items = prns
            pdc = props['sub2cls']
        elif key == 'ObjectProperty':
            items = oppys
            pdc = props['ppy2prn']
        else:
            items = appys
            pdc = props['ppy2prn']
        pns = pdc[0].lower()
        pdcU = '<' + uris[pns] + pdc[1] + '>'
        lpdc = aprops['sth2lbl']
        lns = lpdc[0].lower()
        lpdcU = '<' + uris[lns] + lpdc[1] + '>'
        for item in items:
            root = pydic[key][item]
            ins = root[0].lower()
            itemU = '<' + uris[ins] + root[1] + '>'
            itemUs[item] = itemU
            objU = '<' + uris['rdfs'] + key + '>'
            buff += itemU + ' ' + pdcU + ' ' + objU + END
            nm = DQ + root[2] + DQ
            buff += itemU + ' ' + lpdcU + ' ' + nm + END
    return(buff)

def bridge(itemUs, stmU, clsU, pdcU, lftU, rhtU, mynm=''):
    buff = ''
    rdfBU = uris['rdf']
    cnt = 1
    if mynm:
        buff += clsU + ' ' + itemUs['sub2cls'] + ' ' + stmU + END
        buff += clsU + ' ' + itemUs['sth2lbl'] + ' ' + DQ + mynm + DQ + END
        buff += clsU + ' ' + '<' + rdfBU + 'subject' + '>' + ' ' + lftU + END
        buff += clsU + ' ' + '<' + rdfBU + 'object' + '>' + ' ' + rhtU + END
        buff += clsU + ' ' + '<' + rdfBU + 'predicate' + '>' + ' ' + pdcU + END
        cnt += 5
    buff += lftU + ' ' + pdcU + ' ' + rhtU + END
    out = (buff, cnt)
    return(out)

def me2you(item, key, sbjU, pdcU, objBU):
    buff = ''
    try:
        vals = sorted(item[key].keys())
    except KeyError:
        return(0)
    cnt = len(vals)
    for val in vals:
        if(key == 'hgncg'):
            val = val.split(':')[1]
        myU = '<' + objBU + val + '>'
        buff += sbjU + ' ' + pdcU + ' ' + myU + END
    out = (buff, cnt)
    return(out)

def you2me(item, key, objU, pdcU, sbjBU):
    buff = ''
    try:
        vals = sorted(item[key].keys())
    except KeyError:
        return(0)
    cnt = len(vals)
    for val in vals:
        #abits = val.split('-')
        #val = abits[0] # canonocal accession
        if(key == 'hgncg'):
            val = val.split(':')[1]
        myU = '<' + sbjBU + val + '>'
        buff += myU + ' ' + pdcU + ' ' + objU + END
    out = (buff, cnt)
    return(out)

def me2lbl(item, key, sbjU, pdcU, dlm='_', ind=0):
    buff = ''
    try:
        vals = sorted(item[key].keys())
    except KeyError:
        return(0)
    cnt = 0
    for val in vals:
        bits = val.split(dlm)
        try:
            val = bits[ind]
        except KeyError:
            return(0)
        buff += sbjU + ' ' + pdcU + ' ' + DQ + val + DQ + END
        cnt += 1
    out = (buff, cnt)
    return(out)

def upac2xref(upac, mapA):
    parcs = []
    if upac in mapA.keys():
        parcs = sorted(mapA[upac].keys())
    return parcs

def mapAlert(lst, lstlbl, var1, lbl1):
    mycnt = len(lst)
    if mycnt == 1:
        return()
    if mycnt > 1:
        mycnt = str(mycnt)
        msg = 'WARNING: mapAlert: ' + lbl1 + ':' + var1 + ':' + lstlbl + ':' + mycnt + ':' + str(lst)
    elif mycnt == 0:
        mycnt = str(mycnt)
        msg = 'ERROR: mapAlert: ' + lbl1 + ':' + var1 + ':' + lstlbl + ':' + mycnt + ':' + str(lst)
    print(msg)

def xduo(idmdat, xkey, mapA={}):
    data = {
    }
    mapkeys = sorted(mapA.keys())
    for upac in sorted(idmdat.keys()):
        if mapkeys:
            if upac not in mapkeys:
                continue
        idmacc = idmdat[upac]
        upbac = upac.split('-')[0]
        if upac != upbac:
            idmbac = idmdat[upbac]
        try:
            xkeys = sorted(idmacc[xkey].keys())
        except KeyError:
            if upac != upbac:
                try:
                    xkeys = sorted(idmbac[xkey].keys())
                except KeyError:
                    continue
            else:
                continue
        # multiple gene names (not synonyms) for an accession (basic + iso)
        mapAlert(xkeys, 'xduo.xkeys', upac, 'xduo.upac')
        for onekey in xkeys:
            if onekey not in data.keys():
                data[onekey] = {}
            if upac not in data[onekey].keys():
                data[onekey][upac] = {}
            data[onekey][upac] = idmacc

    return(data)

def upac2bac(idmdat, mapA={}):
    data = {
    }
    mapkeys = sorted(mapA.keys())
    for upac in sorted(idmdat.keys()):
        if mapkeys:
            if upac not in mapkeys:
                continue
        idmacc = idmdat[upac]
        upbac = upac.split('-')[0]
        if upac != upbac:
            idmbac = idmdat[upbac]
        if upbac not in data.keys():
            data[upbac] = {}
        if upac not in data[upbac].keys():
            data[upbac][upac] = {}
        data[upbac][upac] = idmacc

    return(data)
############################# global vars #####################################

# not yet used
#'pur' : ['sio', 'SIO_010295', 'process up-regulation'],
#'pdr' : ['sio', 'SIO_010296', 'process down-regulation'],
#'rgts' : ['sio', 'SIO_001125', 'regulation of transcription'],
#'pcx' : ['obo', 'GO_0032991', 'protein complex'],
# 'mrna' : ['sio', 'SIO_010099', 'messanger RNA'],
# 'cc' : ['obo', 'GO_0005575', 'cellular component'], # > 1 gene product!
# 'mf' : ['obo', 'GO_0003674', 'molecular function'],
# 'gp' : ['NCIt', 'C26548', 'gene product'], # including various combinations
# 'pyp' : ['obo', 'CHEBI', '15841', 'polypeptide'],
# depricated
# 'role' : ['sio', 'SIO_000016', 'role'],
# 'rtn' : ['obo', 'GO_0000122', 'negative regulation of transcription from RNA polymerase II promoter'],
# 'ice' : ['sio', 'SIO_010015', 'information content entity'],
# 'rtp' : ['obo', 'GO_0045944', 'positive regulation of transcription from RNA polymerase II promoter'],
#'prl2prl' : ['sio', 'SIO_000630', 'is paralogous to'],
# added 2016
# added 2018
#'gp2cc' : ['sio', 'SIO_000093', 'is proper part of'], # to be used in GOA
#'gp2bp' : ['sio', 'SIO_000062', 'is participant in'], # to be used in GOA # parent of 'is agent in'
#'tlp2tlp' : ['sio', 'SIO_000203', 'is connected to'], # to be used in Intact
#'mbr2lst' : ['sio', 'SIO_000095', 'is member of'],
#'sth2sth' : ['sio', 'SIO_000001', 'is related to'],
#'sth2sbj' : ['sio', 'SIO_000139', 'has agent'],
#'sth2obj' : ['sio', 'SIO_000291', 'has target'],
## not yet used
#'cc2gp' : ['sio', 'SIO_000053', 'has proper part'],
#'sth2pvd' : ['schema', 'provider', 'has provider'],
## depricated
#'rgr2trg' : ['sio', 'SIO_001154', 'regulates'], # tftg
#'acr2trg' : ['sio', 'SIO_001401', 'positively regulates'], # tftg
#'spr2trg' : ['sio', 'SIO_001402', 'negatively regulates'], # tftg
# 'rgr2trg' : ['obo', 'RO_0002448', 'molecularly controls'], # tftg
# 'acr2trg' : ['obo', 'RO_0002450', 'molecularly increases activity of'], # tftg
# 'spr2trg' : ['obo', 'RO_0002449', 'molecularly decreases activity of'], # tftg
#'gn_up' : 'http://rdf.biogateway.eu/gene-uniprot/',
#'tlp2ptm': ['obo', 'RO_0000053', 'bearer of'], # e.g. protein -> modified residue, the same semantics as BFO # NOT 'all-some'
#'cls2cls': ['owl', 'equivalentClass', 'is equivalent class of'],
#'sth2rlm': ['skos', 'relatedMatch', 'has related match'],
#'ppi2tlp': ['sio', 'SIO_000139', 'has agent'], # Attn: by definition between a process and an entity ! Not in BFO
#'bp2gp': ['sio', 'SIO_000132', 'has participant'], #  parent of 'has agent;
olders = {
'bag': ['rdf', 'Bag', 'unordered collection'],
'stm': ['rdf', 'Statement', 'triple'],
'tlp': ['sio', 'SIO_010043', 'protein'],
'gn': ['sio', 'SIO_010035', 'gene'],
'chr': ['obo', 'GO_0005694', 'chromosome'],
'gom': ['obo', 'SO_0001026', 'genome'],
#'ptm': ['obo', 'PR_000025513', 'modified amino-acid residue'],
#'txn': ['sio', 'SIO_010000', 'organism'],
#'ppi': ['obo', 'INO_0000311', 'protein-protein interaction'],
#'tsf': ['ncit', 'C17207', 'transcription factor'],
#'tsfrx': ['obo', 'GO_0090575', 'RNA polymerase II transcription factor complex'],
#'dss': ['sio', 'SIO_010299', 'disease'],
#'rgts': ['obo', 'GO_0006357', 'regulation of transcription from RNA polymerase II promoter'],
#'bp': ['obo', 'GO_0008150', 'biological process'],
#'mi': ['obo', 'MI_0000', 'molecular interaction'],
}
props = {
## in all:
'sub2cls': ['rdfs', 'subClassOf', 'is subclass of'], # used in all
'ppy2prn': ['rdfs', 'subPropertyOf', 'is subproperty of'], # used in header()
'sth2evd': ['sio', 'SIO_000772', 'has evidence'], # PubMed only: ALL
'sth2ori': ['schema', 'evidenceOrigin', 'has evidence origin'], # DATA sources; ALL; e.g. bgwp 2 upca
'sth2src': ['sio', 'SIO_000253', 'has source', 'has source is a relation between an entity and another entity from which it stems from.'], # TODO sth2src => iaf2src (e.g. graph, statement)
## only in entities
'gn2txn': ['obo', 'BFO_0000052', 'inheres in'], # only in gene.py, TODO expunge
'gn2chr': ['obo', 'BFO_0000050', 'part of'],
'sub2set': ['obo', 'BFO_0000050', 'part of'],
'gp2txn': ['obo', 'BFO_0000052', 'inheres in'], # gene or gene product
'sth2clm' : ['skos', 'closeMatch', 'has close match'], # gene.py grot.py
'mbr2lst': ['schema', 'memberOf', 'is member of'], # gene.py grot.py
# only in bridges
'ins2cls': ['rdf', 'type', 'is instance of'], # used in all bridges
'gp2phn': ['obo', 'RO_0002331', 'involved in'], # not in SIO or BFO2
'gp2bp': ['obo', 'RO_0002331', 'involved in', 'biological_process'], # comes from GPA # TODO fix this quick fix ?
'gp2cc': ['obo', 'BFO_0000050', 'part of', 'cellular_component'], # comes from GPA, present in RO
'gp2mf': ['obo', 'RO_0002327', 'enables', 'molecular_function'], # comes from GPA # not in SIO
'gn2gp': ['sio', 'SIO_010078', 'encodes'], # only in gene.py
'tlp2tlp': ['obo', 'RO_0002436', 'molecularly interacts with'], # only in up2up.py
'rgr2trg': ['obo', 'RO_0002428', 'involved in regulation of'], # only in tsf2gn.py
'sth2mtd': ['rdfs', 'isDefinedBy', 'is defined by'], # up2go.py up2up.py
'orl2orl' : ['sio', 'SIO_000558', 'is orthologous to'],
## TODO to be used
'bml2txn': ['obo', 'IAO_0000115', 'produced by'], # to be used for anything secreted
'acr2trg': ['obo', 'RO_0002429', 'involved in positive regulation of'], # to be used
'spr2trg': ['obo', 'RO_0002430', 'involved in negative regulation of'], # to be used
}
aprops = {
## in all:
'sth2lbl': ['skos', 'prefLabel', 'has name'], # in all
## only in entities
'sth2dfn': ['skos', 'definition', 'has definition'], # in gene.py prot.py
'sth2syn': ['skos', 'altLabel', 'has synonym'], # in gene.py prot.py
'evd2lvl': ['schema', 'evidenceLevel', 'has evidence level'], # prot.py tsf2gn.py up2up.py
# only in bridges
'sth2val': ['rdf', 'value', 'has value'], # only tsf2gn.py
'sth2cmt': ['rdfs', 'comment', 'has comment'], # up2mim.py
## to be used
'sth2id': ['skos', 'notation', 'has notation'], # not used yet
}
pydic = {
'ObjectProperty': props,
'AnnotationProperty': aprops,
'Class': olders,
}
xuris = {
'ncit': 'http://identifiers.org/ncit/',
'upbac': 'http://identifiers.org/uniprot/', # strictly accs
'ncbipr': 'ncbi.nlm.nih.gov/protein/', # accepts both UPID and UPAC
'ncbigi': 'http://identifiers.org/ncbigi/gi:', # the same id for genes and prots, better avoided
'ccds': 'http://identifiers.org/ccds/',
'ncbiug': 'http://www.ncbi.nlm.nih.gov/unigene',
'idungn': 'http://identifiers.org/unigene/', # collection of transcripts associated with a gene
'idomim': 'http://identifiers.org/omim/',
'kegg': 'http://identifiers.org/kegg.genes/',
'embl': 'http://identifiers.org/ena.embl/',
'enspr': 'http://identifiers.org/ensembl.protist/',
'ensfu': 'http://identifiers.org/ensembl.fungi/',
'ensme': 'http://identifiers.org/ensembl.fungi/',
'enspl': 'http://identifiers.org/ensembl.plant/',
'ens': 'http://identifiers.org/ensembl/',
'genewiki': 'http://identifiers.org/genewiki/', # only human
'hgncnm': 'http://identifiers.org/hgnc.symbol/',
'hgncg': 'http://identifiers.org/hgnc/',
'pr': 'http://purl.obolibrary.org/obo/PR_',
'upac': 'http://purl.obolibrary.org/obo/PR_',
}
uris = {
'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
'owl': 'http://www.w3.org/2002/07/owl#',
'skos': 'http://www.w3.org/2004/02/skos/core#',
'schema': 'http://schema.org/',
'obo': 'http://purl.obolibrary.org/obo/',
'sio': 'http://semanticscience.org/resource/',  # resolvable (2014-09-29)
'uniprot': 'http://uniprot.org/uniprot/', # accepts both UPID and UPAC
'uparc': 'http://identifiers.org/uniparc/',
'ncbig': 'http://identifiers.org/ncbigene/',
'ncbitx': 'http://purl.bioontology.org/ontology/NCBITAXON/',
'rfsq': 'http://identifiers.org/refseq/',
'keggortho': 'http://identifiers.org/kegg.orthology/',
#'orthodb': 'https://www.orthodb.org/?query=', # accepts IDs from UP idmapping
'orthodb': 'https://www.orthodb.org/',
'omim': 'http://purl.bioontology.org/ontology/OMIM/',
'intact': 'http://identifiers.org/intact/',
'goa': 'http://identifiers.org/goa/',
'pubmed': 'http://identifiers.org/pubmed/',
'ensp': 'http://identifiers.org/ensembl/',
'ensg': 'http://identifiers.org/ensembl/',
'go': 'http://purl.obolibrary.org/obo/GO_',
'bgw': 'http://rdf.biogateway.eu/',
'gene': 'http://rdf.biogateway.eu/gene/',
'prot': 'http://rdf.biogateway.eu/prot/',
'tsf_gn': 'http://rdf.biogateway.eu/prot-gene/', # for TF-TG
'up_up': 'http://rdf.biogateway.eu/prot-prot/',
'up_obo': 'http://rdf.biogateway.eu/prot-obo/',
'up_mim': 'http://rdf.biogateway.eu/prot-omim/',
'extri': 'http://www.extri.org',
'htri': 'http://www.lbbc.ibb.unesp.br/htri',
'signor': 'http://signor.uniroma2.it',
'tfacts': 'http://www.tfacts.org',
'trrust': 'http://www.grnpedia.org/trrust/',
}
zeno = {
'Opys': props,
'Apys': aprops,
'Prns': olders,
'Uris': uris,
}

with open('../../json/zeno.json', 'w') as outfile:
    json.dump(zeno, outfile, indent=2)

