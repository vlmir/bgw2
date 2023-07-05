#! /usr/bin/env python3
import sys
import glob
import json
import pprint
import zeno as Z
import time

def dat2rdf(stms, bgwmap, path):
    olders = [
    'stm'
    ]
    ### Object properties ###
    props = [
    'gp2phn',
    'sth2ori',
    'sth2src',
    'ins2cls',
    'cls2prn',
    ]
    ### Annotation properties ###
    aprops = [
    'sth2cmt',
    'sth2nm',
    ]
    itemUs = {'self': '<http://rdf.biogateway.eu/graph/prot2phen>'}
    buff = Z.header(itemUs, props, aprops, olders)
    fp = open(outpth, 'w')
    fp.write(buff)
    buff = ''
    ###########################################################################

    protBU = Z.uris['prot']
    omimBU = Z.uris['omim']
    pdcU = itemUs['gp2phn']
    ns = 'up_mim'
    myBU = Z.uris[ns] 
    
    mysrc = 'uniprot'
    srcBU = Z.uris[mysrc] # do NOT delete!
    srcU = '<' + srcBU + '>'
    buff += itemUs['self'] + ' ' + itemUs['sth2src'] + ' ' + srcU + Z.END
    outcnt = 1
    stmU = itemUs['stm']
    txn = '9606'
    for stmk in sorted(stms.keys()):
        stm = stms[stmk]
        (lftid, rhtid) = stmk.split(Z.KS)
        upacL = lftid
        try:
            ourids = sorted(bgwmap['upac'][upacL]['bgw'].keys())
        except KeyError:
            print("KeyError:dat2rdf:bgwmap['upac']:" + upacL)
            continue
        if len(ourids) != 1:
            Z.mapAlert(ourids, 'ourids', upacL, 'upac')
        rhtU = '<' + omimBU + rhtid + '>'
        for lftid in ourids:
            lftU = '<' + protBU + lftid + '>'
            ###
            mynm = 'protein - disease interaction ' + upacL + ':' + rhtid
            lftid = lftid.replace('/', Z.SS) # Attn: renaming
            clsid = lftid + Z.US + rhtid
            clsU = '<' + myBU + clsid + '>'
            insid = lftid + Z.US + rhtid + '#' + mysrc
            insU = '<' + myBU + insid + '>'
            out = Z.bridge(itemUs, stmU, clsU, pdcU, lftU, rhtU, mynm)
            if out:
                buff += out[0]
                outcnt += out[1]
            outcnt += 6
            buff += insU + ' ' + itemUs['ins2cls'] + ' ' + clsU + Z.END
            outcnt += 1
            #buff += insU + ' ' + itemUs['sth2ori'] + ' ' + srcU + Z.END
            ## source
            onesrcU = '<' + srcBU + upacL + '>'
            mypdcU = itemUs['sth2ori']
            buff += insU + ' ' + mypdcU + ' ' + onesrcU + Z.END
            outcnt += 1
    
            ###
            if 'change' in stm.keys(): # TODO ALWAYS use this construct!
                for key in sorted(stm['change'].keys()):
                    cmt = Z.DQ + key + Z.DQ
                    buff += insU + ' ' + itemUs['sth2cmt'] + ' ' + cmt + Z.END
                    outcnt += 1

        fp.write(buff) # per upac
        buff = ''

    fp.close()
    print(outpth + ':' + str(outcnt))

    return(0) 
    ########################### dat2rdf() ####################################

## AARS      P49588     VAR_063527  p.Arg329His    Disease       rs267606621 Charcot-Marie-Tooth disease 2N(CMT2N) [MIM:613287]
##  1 gene name
## 11 up acc
## 22
## 34 aa change
## 49 type
## 63 '-' possible
## 75 description

txtkeys = {
'uniprot':[10, 20],
'change':[33, 47],
'omim':[74, -1]
}

dirpth = sys.argv[1] 
txt = sys.argv[2]
txts = glob.glob(dirpth + '*.' + txt)
xrf = sys.argv[3]
print(time.ctime(time.time()))

for txtpth in txts:
    txtdat = {}
    txtdat = Z.xfixlen(txtpth, txtkeys)
    count = len(txtdat.keys())
    if not count:
        print('WARNING: no data in: ' + txtpth)
        continue
    print(str(txtpth) + ':' + str(count))
    xrfpth = txtpth.replace(txt, xrf)
    with open(xrfpth, "r") as fp:
        bgwmap = json.load(fp)

    #pprint.pprint(txtdat)
    outpth = txtpth.replace(txt, 'nt')
    dat2rdf(txtdat, bgwmap, outpth)
print(time.ctime(time.time()))
