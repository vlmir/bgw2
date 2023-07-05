#! /usr/bin/env python3
import sys
import glob
import json
import pprint
import csv
import zeno as Z
import time

def dat2rdf(alldat, bgwmap, gxrfdat, outpth):
    olders = [
    'stm'
    ]
    ### Object properties ###
    props = [
    #'acr2trg',
    #'spr2trg',
    'rgr2trg',
    'sth2ori',
    'sth2evd',
    'ins2cls',
    'cls2prn',
    ]
    ### Annotation properties ###
    aprops = [
    'evd2lvl',
    'sth2val',
    'sth2nm',
    ]
    itemUs = {}
    itemUs = {'self': '<http://rdf.biogateway.eu/graph/tfac2gene>'}
    buff = Z.header(itemUs, props, aprops, olders)
    fp = open(outpth, 'w')
    fp.write(buff)
    buff = ''
    ###########################################################################

    pdcU = itemUs['rgr2trg'] # occurs on TFacts
    clsU = itemUs['stm']
    ns = 'tsf_gn'
    myBU = Z.uris[ns] 
    protBU = Z.uris['prot']
    geneBU = Z.uris['gene']
    stmU = itemUs['stm']
    outcnt = 0
    bdat = alldat['base']
    for key in sorted(bdat.keys()): # key: tftg id
        kbdat = bdat[key] # base/common data
        for lftid in sorted(kbdat['accs'].keys()):
            upacL = lftid
            try:
                ourids = sorted(bgwmap['upac'][upacL]['bgw'].keys())
            except KeyError: # ucount 3
                # the line below commented for testing with p53 alone
                #print("KeyError:dat2rdf:bgwmap['upac'][" + upacL + "]['bgw']")
                continue
            if len(ourids) != 1:
                Z.mapAlert(ourids, 'ourids', upacL, 'upac')
            for lftid in ourids:
                lftU = '<' + protBU + lftid + '>'
                lftid = lftid.replace('/', Z.SS) # Attn: renaming

                for ncbig in sorted(kbdat['ncbig'].keys()):
                    try:
                        ourids = sorted(gxrfdat['ncbig'][ncbig]['bgw'].keys())
                    except KeyError: # ucount:624, seem absent in the human UP files
                        #print("KeyError:dat2rdf:gxrfdat['ncbig']:" + ncbig)
                        continue
                    for rhtid in ourids:
                        rhtU = '<' + geneBU + rhtid + '>'
                        ### exporting
                        ## basic data
                        rhtid = rhtid.replace('/', Z.SS) # Attn: renaming
                        rhtid = rhtid[:-1] # removing the last character
                        clsid = lftid + Z.US + rhtid
                        clsU = '<' + myBU + clsid + '>'
                        mynm = 'transcription factor - target gene interaction ' + key
                        out = Z.bridge(itemUs, stmU, clsU, pdcU, lftU, rhtU, mynm)
                        if out:
                            buff += out[0]
                            outcnt += out[1]
                        ## instances
                        for src in sorted(srckeys.keys()):
                            srclc = src.lower()
                            srcU = '<' + Z.uris[srclc] + '>'
                            try:
                                ksdat = alldat[src][key] # source specific
                            except KeyError:
                                continue
                            if not ksdat['srcs']: # SIC!! - real check for data
                                continue
                            insid = lftid + Z.US + rhtid + '#' + srclc
                            insU = '<' + myBU + insid + '>'
                            buff += insU + ' ' + itemUs['ins2cls'] + ' ' + clsU + Z.END
                            buff += insU + ' ' + itemUs['sth2ori'] + ' ' + srcU + Z.END
                            outcnt += 2
                            # publications
                            mykey = 'refs'
                            if mykey in ksdat.keys():
                                for chunk in ksdat['refs'].keys():
                                    pdcUk = itemUs['sth2evd']
                                    if chunk == '':
                                        continue
                                    if chunk == ' ':
                                        continue
                                    try:
                                        int(chunk)
                                    except ValueError: # occurs on TFacts
                                        continue
                                    refU = '<' + Z.uris['pubmed'] + chunk + '>'
                                    buff += insU + ' ' + pdcUk + ' ' + refU + Z.END
                                    outcnt += 1
                            # conficence level
                            mykey = 'cnfs'
                            if mykey in ksdat.keys():
                                for chunk in ksdat['cnfs'].keys():
                                    pdcUk = itemUs['evd2lvl']
                                    vals = []
                                    maxval = 0
                                    # scaling factor (2i+1)/2N*100
                                    # N:2, i in (0,N-1)
                                    if chunk == 'Low':
                                        chunk = 25
                                    elif chunk == 'High':
                                        chunk = 75
                                    else:
                                        chunk = 0
                                    vals.append(chunk)
                                    maxval = max(vals)
                                    obj = Z.DQ + str(maxval) + Z.DQ
                                    buff += insU + ' ' + pdcUk + ' ' + obj + Z.END
                                    outcnt += 1
                            # regulation mode
                            mykey = 'signs'
                            if mykey in ksdat.keys():
                                for chunk in ksdat[mykey].keys():
                                    if(chunk == 'UP'):
                                        mode = Z.DQ + 'positive' + Z.DQ
                                        buff += insU + ' ' + itemUs['sth2val'] + ' ' + mode + Z.END
                                        outcnt += 1
                                    if(chunk == 'DOWN'):
                                        mode = Z.DQ + 'negative' + Z.DQ
                                        buff += insU + ' ' + itemUs['sth2val'] + ' ' + mode + Z.END
                                        outcnt += 1

        fp.write(buff) # per gene
        buff = ''

    fp.close()
    print(outpth + ':' + str(outcnt))

    return(0)

######################### END OF dat2rdf() ###################################

bkeys = {
0:'fgid',
1:'accs',
2:'ncbig',
32:'tfnm',
33:'tgnm',
37:'tgensg',
}
srckeys = {
'EXTRI':{
5:'srcs',
4:'refs',
3:'cnfs'
},
'HTRI':{
6:'srcs',
8:'refs',
9:'cnfs'
},
'TRRUST':{
10:'srcs',
11:'signs',
12:'refs',
},
'TFacts':{
13:'srcs',
14:'signs',
17:'refs',
18:'cnfs'
},
'GOA':{
19:'srcs',
20:'signs',
},
'Intact':{
21:'srcs',
22:'refs',
},
'Signor':{
24:'srcs',
26:'signs',
27:'refs',
}
}
dirpth = sys.argv[1]
ori = sys.argv[2]
xrf = sys.argv[3]
gxrf = sys.argv[4]
oris = glob.glob(dirpth + '*.' + ori)
print(time.ctime(time.time()))
for oripth in oris:
    alldat = {
    'base':{}
    }
    bdat = Z.x2D(oripth, '\t', '|', bkeys, 0)
    count = len(bdat.keys())
    if not count:
        print('WARNING: no data in: ' + oripth)
        continue
    print(str(oripth) + ':' + str(count))
    xrfpth = oripth.replace(ori, xrf)
    with open(xrfpth, "r") as fp:
        bgwmap = json.load(fp)
    gxrfpth = oripth.replace(ori, gxrf)
    with open(gxrfpth, "r") as fp:
        gxrfdat = json.load(fp)
    alldat['base'] = bdat
    outpth = oripth.replace(ori, 'nt')
    for src in sorted(srckeys.keys()):
        skeys = srckeys[src] # dict
        sdat = Z.x2D(oripth, '\t', ';', skeys, 0)
        count = len(sdat.keys())
        if not count:
            continue
        if src not in alldat.keys():
            alldat[src] = {}
        alldat[src] = sdat
    dat2rdf(alldat, bgwmap, gxrfdat, outpth)
print(time.ctime(time.time()))
