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
    refs = sorted(data['refs'].keys())
    for ref in refs:
        (rdb, rid) = ref.split(':')
        if rdb != 'PMID':
            continue
        refU = '<' + Z.uris['pubmed'] + rid + '>'
        buff += insU + ' ' + itemUs['sth2evd'] + ' ' + refU + Z.END
        outcnt += 1
    ecos = sorted(data['eco'].keys())
    ecocnt = len(ecos)
    if ecocnt < 1:
        # nimerous multiple ECO IDs
        print('WARNING: meta: unexpected number of ECO IDs: ' + str(ecocnt) + ' for: ' + insU)
    for eco in ecos:
        ecoid = eco.replace(':', '_')
        ecoU = '<' + Z.uris['obo'] + ecoid + '>'
        buff += insU + ' ' + itemUs['sth2mtd'] + ' ' + ecoU + Z.END
        outcnt += 1
    outlst = [buff, outcnt]
    return(outlst)

def dat2rdf(gpadat, bgwmap, outpth):
    olders = [
    'stm',
    ]
    ### Object properties ###
    props = [
    'cls2prn',
    'gp2bp',
    'gp2cc',
    'gp2mf',
    'ins2cls',
    'sth2evd',
    'sth2mtd',
    'sth2ori',
    'sth2src',
    ]
    ### Annotation properties ###
    aprops = [
    'sth2nm',
    ]
    itemUs = {'self': '<http://rdf.biogateway.eu/graph/prot2onto>'}
    buff = Z.header(itemUs, props, aprops, olders)
    fp = open(outpth, 'w')
    fp.write(buff)
    buff = ''
    ###########################################################################

    rdfU = Z.uris['rdf']
    pdcks = {
    'part_of':'gp2cc',
    'enables':'gp2mf',
    'involved_in':'gp2bp',
    }
    gosubs = {
    'part_of':'cellular component',
    'enables':'molecular finction',
    'involved_in':'biological process',
    }
    clsU = itemUs['stm']
    ns = 'up_obo'
    myBU = Z.uris[ns]
    protBU = Z.uris['prot']
    oboBU = Z.uris['obo']
    stmU = itemUs['stm']
    mysrc = 'goa'
    srcBU = Z.uris[mysrc] # do NOT delete!
    srcU = '<' + srcBU + '>'
    buff += itemUs['self'] + ' ' + itemUs['sth2src'] + ' ' + srcU + Z.END
    outcnt = 1
    for lftid in sorted(gpadat.keys()): # UPAC
        lftdat = gpadat[lftid]
        #lftid = lftid.replace('PRO_', '') # PR bare ID TODO implement
        shreds = lftid.split(':')
        if len(shreds) == 2: # AC:PRO_
            (upacL, ext) = shreds
        else:
            upacL = lftid
        try:
            ourids = sorted(bgwmap['upac'][upacL]['bgw'].keys())
        except KeyError: # count:1243 TODO
            #print("KeyError:dat2rdf:bgwmap['upac']:" + upacL)
            continue
        if len(ourids) != 1:
            Z.mapAlert(ourids, 'ourids', upacL, 'upac')
        for lftid in ourids:
            lftU = '<' + protBU + lftid + '>'
            for rhtid in sorted(lftdat.keys()): # the value '-' occurs for some taxa
                rhtdat = lftdat[rhtid]
                rhtid = rhtid.replace(':', '_')
                rhtbid = rhtid
                qlrs = sorted(rhtdat['qlrs'].keys())
                myqlrs = {}
                for qlr in qlrs:
                    if qlr in pdcks.keys(): # 'colocalizes_with' occurs
                        pdck = pdcks[qlr]
                        myqlrs[qlr] = 1
                    else:
                        continue
                if len(myqlrs.keys()) != 1:
                    continue
                myqlr = sorted(myqlrs.keys())[0]
                pdcU = itemUs[pdck]
                rhtU = '<' + oboBU + rhtid + '>'
                ### exporting
                ## basic data
                lftid = lftid.replace('/', Z.SS) # Attn: renaming
                clsid = lftid + Z.US + rhtid
                clsU = '<' + myBU + clsid + '>'
                mynm = 'protein - ' + gosubs[myqlr] + ' association ' + upacL + ':' + rhtid
                out = Z.bridge(itemUs, stmU, clsU, pdcU, lftU, rhtU, mynm)
                if out:
                    buff += out[0]
                    outcnt += out[1]
                insid = lftid + Z.US + rhtid + '#' + mysrc
                insU = '<' + myBU + insid + '>'
                buff += insU + ' ' + itemUs['ins2cls'] + ' ' + clsU + Z.END
                #buff += Z.triple(insU, lftU, rhtU, pdcU)
                buff += insU + ' ' + itemUs['sth2ori'] + ' ' + srcU + Z.END
                outcnt += 2
                ## metadata
                outlst = meta(insU, rhtdat, itemUs)
                buff += outlst[0]
                outcnt += outlst[1]

        fp.write(buff) # per upac
        buff = ''

    fp.close()
    print(outpth + ':' + str(outcnt))

    return(0)
################################ dat2rdf() ####################################

gpakeys = {
0:'db',
1:'acc',
2:'qlrs',
3:'goid',
4:'refs', # not all have pubmed refs!
5:'eco', # single
}

dirpth = sys.argv[1] 
gpa = sys.argv[2]
gpas = glob.glob(dirpth + '*.' + gpa)
xrf = sys.argv[3]
print(time.ctime(time.time()))

for gpapth in gpas:
    # Note: '|' seems to be used in gpa files only in the 7th field
    gpadat = Z.x2Dx2D(gpapth, '\t', '|', gpakeys, 1, 3) # no filtering
    count = len(gpadat.keys())
    if not count:
        print('WARNING: no gpadat in: ' + gpapth)
        continue
    print(str(gpapth) + ':' + str(count))
    xrfpth = gpapth.replace(gpa, xrf)
    with open(xrfpth, "r") as fp:
        bgwmap = json.load(fp)

    #pprint.pprint(gpadat)
    outpth = gpapth.replace(gpa, 'nt')
    dat2rdf(gpadat, bgwmap, outpth)
print(time.ctime(time.time()))
