#! /usr/bin/env python3
import sys
import glob
import json
import pprint
import zeno as Z
import time

def meta(insU, data, itemUs):
    buff = ''
    outcnt = 0
    for key in sorted(data['pubids'].keys()):
        bits = key.split(':') # very many refs 'psi-mi:MI:...'
        if bits[0] == 'pubmed':
            kid = bits[1]
            itemU = '<' + Z.uris['pubmed'] + kid + '>'
            buff += insU + ' ' + itemUs['sth2evd'] + ' ' + itemU + Z.END
            outcnt += 1
    for key in sorted(data['cnfvals'].keys()):
        bits = key.split(':') # very many 'author score:...'
        if bits[0] == 'intact-miscore':
            kid = str(int(100*float(bits[1])))
            buff += insU + ' ' + itemUs['evd2lvl'] + ' ' + Z.DQ + kid + Z.DQ + Z.END
            outcnt += 1
    for key in sorted(data['inacids'].keys()):
        bits = key.split(':') # very many 'psi-mi:MI:...'
        if bits[0] == 'intact':
            kid = bits[1]
            itemU = '<' + Z.uris['intact'] + kid + '>'
            buff += insU + ' ' + itemUs['sth2ori'] + ' ' + itemU + Z.END
            outcnt += 1
    for key in sorted(data['inactps'].keys()):
        bits = key.split(Z.DQ)
        if bits[0] == 'psi-mi:':
            kid = bits[1].replace(':', '_')
            itemU = '<' + Z.uris['obo'] + kid + '>'
            buff += insU + ' ' + itemUs['ins2cls'] + ' ' + itemU + Z.END
            outcnt += 1
    for key in sorted(data['mtds'].keys()):
        bits = key.split(Z.DQ)
        if bits[0] == 'psi-mi:':
            kid = bits[1].replace(':', '_')
            itemU = '<' + Z.uris['obo'] + kid + '>'
            buff += insU + ' ' + itemUs['sth2mtd'] + ' ' + itemU + Z.END
            outcnt += 1
    outlst = [buff, outcnt]
    return(outlst)

def dat2rdf(tsvdat, bgwmap, outpth):
    olders = [
    'stm',
    ]
    ### Object properties ###
    props = [
    'cls2prn',
    'ins2cls',
    'sth2evd',
    'sth2mtd',
    'sth2ori',
    'sth2src',
    'tlp2tlp',
    ]
    ### Annotation properties ###
    aprops = [
    'evd2lvl',
    'sth2nm',
    ]
    itemUs = {'self': '<http://rdf.biogateway.eu/graph/prot2prot>'}
    buff = Z.header(itemUs, props, aprops, olders)
    fp = open(outpth, 'w')
    fp.write(buff)
    buff = ''
    ###########################################################################

    protBU = Z.uris['prot']
    pdcU = itemUs['tlp2tlp']
    ns = 'up_up'
    myBU = Z.uris[ns] 
    stmU = itemUs['stm']
    mysrc = 'intact' 
    srcBU = Z.uris[mysrc] # do NOT delete!
    srcU = '<' + srcBU + '>'
    buff += itemUs['self'] + ' ' + itemUs['sth2src'] + ' ' + srcU + Z.END
    outcnt = 1
    for lftk in sorted(tsvdat.keys()):
        ## filtering
        bits = lftk.split(':')
        if len(bits) != 2:
            continue
        dat1 = tsvdat[lftk]
        (lftdb, lftid) = (bits[0], bits[1])
        if lftdb != 'uniprotkb':
            continue
        #lftid = lftid.replace('PRO_', '') # PR bare ID TODO implement
        shreds = lftid.split(':')
        if len(shreds) == 2: # AC:PRO_
            (upacL, ext) = shreds
        else:
            upacL = lftid
        try:
            ourids = sorted(bgwmap['upac'][upacL]['bgw'].keys())
        except KeyError: # count:63970
            #print("KeyError:dat2rdf:bgwmap['upac']:" + upacL)
            continue
        if len(ourids) != 1:
            Z.mapAlert(ourids, 'ourids', upacL, 'upac')
        for lftid in ourids:
            lftU = '<' + protBU + lftid + '>'
            for rhtk in sorted(dat1.keys()): # the value '-' occurs for some taxa
                bits = rhtk.split(':')
                if len(bits) != 2:
                    continue
                dat2 = dat1[rhtk]
                (rhtdb, rhtid) = (bits[0], bits[1])
                if rhtdb != 'uniprotkb':
                    continue
                shreds = rhtid.split(':')
                if len(shreds) == 2: # AC:PRO_
                    (upacR, ext) = shreds
                else:
                    upacR = rhtid
                try:
                    ourids = sorted(bgwmap['upac'][upacR]['bgw'].keys())
                except KeyError:
                    #print("KeyError:dat2rdf:bgwmap['upac']:" + upacR)
                    continue
                if len(ourids) != 1:
                    Z.mapAlert(ourids, 'ourids', upacR, 'upac')
                for rhtid in ourids:
                    if rhtid < lftid: # sorting ids:
                        swp = lftid
                        lftid = rhtid
                        rhtid = swp
                        swp = lftid
                        lftid = rhtid
                        rhtid = swp
                    rhtU = '<' + protBU + rhtid + '>'
                    ### exporting
                    lftid = lftid.replace('/', Z.SS) # Attn: renaming
                    rhtid = rhtid.replace('/', Z.SS) # Attn: renaming
                    clsid = lftid + Z.US + rhtid
                    clsU = '<' + myBU + clsid + '>'
                    mynm = 'protein - protein interaction ' + upacL + ':' + upacR
                    out = Z.bridge(itemUs, stmU, clsU, pdcU, lftU, rhtU, mynm)
                    if out:
                        buff += out[0]
                        outcnt += out[1]
                    out = Z.bridge(itemUs, stmU, clsU, pdcU, rhtU, lftU) # inverting triple
                    if out:
                        buff += out[0]
                        outcnt += out[1]
                    insid = lftid + Z.US + rhtid + '#' + mysrc
                    insU = '<' + myBU + insid + '>'
                    buff += insU + ' ' + itemUs['ins2cls'] + ' ' + clsU + Z.END
                    #buff += insU + ' ' + itemUs['sth2ori'] + ' ' + srcU + Z.END
                    outcnt += 2
                    ## metadata
                    outlst = meta(insU, dat2, itemUs)
                    buff += outlst[0]
                    outcnt += outlst[1]

        fp.write(buff) # per upac
        buff = ''

    fp.close()
    print(outpth + ':' + str(outcnt))

    return(0)
############################## dat2rdf() #####################################

#     Input file fields:
#            1 ID(s) interactor A
#            2 ID(s) interactor B
#            3 Alt. ID(s) interactor A
#            4 Alt. ID(s) interactor B
#            5 Alias(es) interactor A
#            6 Alias(es) interactor B
#            7 Interaction detection method(s)
#            8 Publication 1st author(s)
#            9 Publication Identifier(s)
#           10 Taxid interactor A
#           11 Taxid interactor B
#           12 Interaction type(s)
#           13 Source database(s)
#           14 Interaction identifier(s)
#           15 Confidence value(s)
#           16 Expansion method(s)
#           17 Biological role(s) interactor A
#           18 Biological role(s) interactor B
#           19 Experimental role(s) interactor A
#           20 Experimental role(s) interactor B
#           21 Type(s) interactor A
#       22 Type(s) interactor B
#       23 Xref(s) interactor A
#       24 Xref(s) interactor B
#       25 Interaction Xref(s)
#       26 Annotation(s) interactor A
#       27 Annotation(s) interactor B
#       28 Interaction annotation(s)
#       29 Host organism(s)
#       30 Interaction parameter(s)
#       31 Creation date
#       32 Update date
#       33 Checksum(s) interactor A
#       34 Checksum(s) interactor B
#       35 Interaction Checksum(s)Negative
#       36 Feature(s) interactor A
#       37 Feature(s) interactor B
#       38 Stoichiometry(s) interactor A
#       39 Stoichiometry(s) interactor B
#       40 Identification method participant A
#       41 Identification method participant B

tsvkeys = {
#0:'idas',
#1:'idbs',
6:'mtds',
8:'pubids',
9:'txnA',
10:'txnB',
11:'inactps',
#12:'srcdbs',
13:'inacids',
14:'cnfvals',
#15:'xpnmtds',
18:'exprlAs',
19:'exprlBs',
28:'hosts',
}

dirpth = sys.argv[1] 
tsv = sys.argv[2]
tsvs = glob.glob(dirpth + '*.' + tsv)
xrf = sys.argv[3]
print(time.ctime(time.time()))

for tsvpth in tsvs:
    tsvdat = Z.x2Dx2D(tsvpth, '\t', '|', tsvkeys, 0, 1) # no filtering
    count = len(tsvdat.keys())
    if not count:
        print('WARNING: no tsvdat in: ' + tsvpth)
        continue
    print(str(tsvpth) + ':' + str(count))
    xrfpth = tsvpth.replace(tsv, xrf)
    with open(xrfpth, "r") as fp:
        bgwmap = json.load(fp)

    #pprint.pprint(tsvdat)
    outpth = tsvpth.replace(tsv, 'nt')
    dat2rdf(tsvdat, bgwmap, outpth)
print(time.ctime(time.time()))
