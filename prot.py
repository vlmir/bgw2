#! /usr/bin/env python3
import sys
import glob
import json
import pprint
import zeno as Z
import time

def dat2rdf(idmdat, tsvdat, outpth):
    olders = [
    'bag',
    'tlp',
    ]
    ### Object properties ###
    props = [
    'gp2txn',
    'sth2ori',
    'sth2evd',
    'sth2src',
    'sth2clm',
    'mbr2lst',
    'cls2prn'
    ]
    ### Annotation properties ###
    aprops = [
    'evd2lvl',
    'sth2nm',
    'sth2syn',
    'sth2dfn',
    ]
    itemUs = {'self': '<http://rdf.biogateway.eu/graph/prot>'}
    buff = Z.header(itemUs, props, aprops, olders)
    fp = open(outpth, 'w')
    fp.write(buff)
    buff = ''
    ###########################################################################
    ### cardinality accession:GeneNmae is N:M !!!

    prnU = itemUs['tlp']
    bagU = '<' + Z.uris['rdf'] + 'Bag' + '>'
    txnBU = Z.uris['ncbitx']
    protBU = Z.uris['prot']
    #prBU = Z.uris['pr'] # do NOT delete!
    mysrc = 'uniprot'
    srcBU = Z.uris[mysrc] # do NOT delete!
    srcU = '<' + srcBU + '>'
    outcnt = 0
    # datafields to be treated as xrfs
    insts = (
    'ensp',
    'rfsq',
    )
    onenmsp = 'upac'
    bgwmap = {
    'bgw': {},
    'upac': {}
    }

    buff += itemUs['self'] + ' ' + itemUs['sth2src'] + ' ' + srcU + Z.END
    idminv = Z.upac2bac(idmdat)
    for onebac in sorted(tsvdat.keys()):
        try:
            idmbac = idminv[onebac]
        except KeyError:
            continue
        idmbac = idminv[onebac]
        unidat = tsvdat[onebac] # Attn: tsvdat going to be modified !!! ???
        unidat[onenmsp] = {} # TODO unidat => bacdat
        upacs = sorted(idmbac.keys())
        # propagating data to isoforms
        for upac in upacs:
            if upac not in unidat['upac'].keys():
                unidat['upac'][upac] = {}
            idmac = idmbac[upac]
            for mynmsp in idmac.keys():
                myids = idmac[mynmsp]
                if myids:
                    unidat['upac'][upac][mynmsp] = {}
                for myid in myids:
                    unidat['upac'][upac][mynmsp][myid] = 1

        ### sanity checks and filtering
        txns = sorted(unidat['txn'].keys())
        if len(txns) != 1: # couunt:0
            Z.mapAlert(txns, 'txns', onebac, onenmsp)
            continue
        gnms = sorted(unidat['gnm'].keys())
        if len(gnms) == 0: # couunt:845, no other alerts
            continue
        if len(gnms) > 1: # couunt:845, no other alerts
            Z.mapAlert(gnms, 'gnms', onebac, onenmsp)
        onegnm = gnms[0]
        mynmsp = 'spnm'
        if mynmsp in unidat.keys():
            spnms = sorted(unidat[mynmsp].keys())
        if len(spnms) == 0: # sic - for some no
            continue
        if len(spnms) > 1: # sic - for some no
            Z.mapAlert(spnms, 'spnms', onegnm, onenmsp)
        mynmsp = 'pome'
        if mynmsp in unidat.keys():
            pomes = sorted(unidat[mynmsp].keys())
        if len(pomes) != 1: # count:99, all muliple, no other alerts
            Z.mapAlert(pomes, 'pomes', onegnm, onenmsp)
        somes = {}
        for pome in pomes:
            bits = pome.split(' ')
            if len(bits) != 3:
                continue # skipping 'Unplaced'
            mypome = bits[0]
            some = bits[2]
            if mypome not in somes.keys():
                somes[mypome] = {}
            somes[mypome][some] = 1
        allpomes = sorted(somes.keys())
        if len(allpomes) == 0:
            continue
        if len(allpomes) > 1:
            Z.mapAlert(allpomes, 'allpomes', onegnm, onenmsp)
        onesome = somes[allpomes[0]]

        ### exporting
        ## defining properties
        txn = txns[0]
        for mysome in onesome.keys():
            if len(mysome) == 1:
                smnm = 'chr-0' + mysome
            else:
                smnm = 'chr-' + mysome
            bagid = txn + '/' + onegnm + '/' # terminating with '/' for bags
            bagid = txn + '/' + smnm + '/' + onegnm + '/' # terminating with '/' for bags
            onebagU = '<' + protBU + bagid + '>'
            mypdcU = itemUs['cls2prn']
            buff += onebagU + ' ' + mypdcU + ' ' + bagU + Z.END
            for upac in upacs:
                mynmsp = 'uparc'
                uparcs = sorted(unidat['upac'][upac][mynmsp].keys())
                if len(uparcs) != 1: # couunt:0
                    Z.mapAlert(uparcs, 'uparcs', upac, onenmsp)
                    continue
                gpid = bagid + uparcs[0]
                gpU = '<' + Z.uris['prot'] + gpid + '>'
                mypdcU = itemUs['cls2prn'] # sic! MUST be re-assigned !!
                buff += gpU + ' ' + mypdcU + ' ' + prnU + Z.END
                outcnt += 1
                # taxon
                mynmsp = 'txn'
                mypdcU = itemUs['gp2txn']
                myBU = txnBU
                out = Z.me2you(unidat, mynmsp, gpU, mypdcU, myBU)
                buff += out[0]
                outcnt += out[1]
    
                ## source
                onesrcU = '<' + srcBU + onebac + '>'
                mypdcU = itemUs['sth2ori']
                buff += gpU + ' ' + mypdcU + ' ' + onesrcU + Z.END
                outcnt += 1
    
                ## name
                bits = upac.split('-')
                ext = '' # needed !!
                if len(bits) > 1:
                    ext = bits[1]
                upids = sorted(unidat['upid'].keys())
                upid = upids[0]
                onenm = upid.split('_')[0]
                if ext:
                    onenm += '-' + ext
                mypdcU = itemUs['sth2nm']
                buff += gpU + ' ' + mypdcU + ' ' + Z.DQ + onenm + Z.DQ + Z.END
                outcnt += 1
    
                ## synonyms
                mypdcU = itemUs['sth2syn']
                # AC
                if onenm != upac:
                    buff += gpU + ' ' + mypdcU + ' ' + Z.DQ + upac + Z.DQ + Z.END
                    outcnt += 1
                # gene name
                if onenm != onegnm:
                    buff += gpU + ' ' + mypdcU + ' ' + Z.DQ + onegnm + Z.DQ + Z.END
                    outcnt += 1
                # gene synonyms
                mynmsp = 'gsnm'
                out = Z.me2lbl(unidat, mynmsp, gpU, mypdcU)
                if out:
                    buff += out[0]
                    outcnt += out[1]
    
                ## description
                mynmsp = 'def'
                if mynmsp in unidat.keys():
                    defs = sorted(unidat[mynmsp].keys())
                if len(defs) != 1: # count:0
                    Z.mapAlert(defs, 'defs', onebac, onenmsp)
                onedef = Z.DS.join(defs) # DS absent in in the output file
                if ext:
                    onedef += ' isoform-' + ext
                mypdcU = itemUs['sth2dfn']
                buff += gpU + ' ' + mypdcU + ' ' + Z.DQ + onedef + Z.DQ + Z.END
                outcnt += 1
    
                ## publications
                mypdcU = itemUs['sth2evd']
                mynmsp = 'pubmed'
                myBU = Z.uris[mynmsp]
                if mynmsp in unidat.keys():
                    out = Z.me2you(unidat, mynmsp, gpU, mypdcU, myBU)
                if out:
                    buff += out[0]
                    outcnt += out[1]
    
                ## annotation scores
                mypdcU = itemUs['evd2lvl']
                mynmsp = 'score'
                if mynmsp in unidat.keys():
                    scores = sorted(unidat[mynmsp].keys())
                    for score in scores:
                        bits = score.split(' ')
                        val = bits[0]
                        buff += gpU + ' ' + mypdcU + ' ' + Z.DQ + val + Z.DQ + Z.END
                        outcnt += 1
    
                ## xrfs
                mypdcU = itemUs['sth2clm']
                myxrfs = unidat['upac'][upac]
                if gpid not in bgwmap['bgw'].keys():
                    bgwmap['bgw'][gpid] = {}
                if upac not in bgwmap['upac'].keys():
                    bgwmap['upac'][upac] = {}
                if 'bgw' not in bgwmap['upac'][upac].keys():
                    bgwmap['upac'][upac]['bgw'] = {}
                if 'upac' not in bgwmap['bgw'][gpid].keys():
                    bgwmap['bgw'][gpid]['upac'] = {}
                bgwmap['upac'][upac]['bgw'][gpid] = 1
                bgwmap['bgw'][gpid]['upac'][upac] = 1
                for mynmsp in insts:
                    if mynmsp not in myxrfs.keys():
                        continue
                    myBU = Z.uris[mynmsp]
                    out = Z.me2you(myxrfs, mynmsp, gpU, mypdcU, myBU)
                    if not out:
                        continue
                    buff += out[0]
                    outcnt += out[1]
                    # building ID map
                    if mynmsp not in bgwmap.keys():
                        bgwmap[mynmsp] = {}
                    if mynmsp not in bgwmap['bgw'][gpid] .keys():
                        bgwmap['bgw'][gpid][mynmsp] = {}
                    for myid in myxrfs[mynmsp]:
                        if myid not in bgwmap[mynmsp].keys():
                            bgwmap[mynmsp][myid] = {}
                            bgwmap[mynmsp][myid]['bgw'] = {}
                        bgwmap[mynmsp][myid]['bgw'][gpid] = 1
                        bgwmap['bgw'][gpid][mynmsp][myid] = 1
    
                # members
                mypdcU = itemUs['mbr2lst']
                buff += gpU + ' ' + mypdcU + ' ' + onebagU + Z.END
                outcnt += 1
    
            fp.write(buff) # per UP canonical AC ?? TODO check it
            buff = ''

    fp.close()
    print(outpth + ':' + str(outcnt))

    return(bgwmap)
    ###################### END OF dat2rdf #####################################

### Fields
idmkeys = {
#'Ensembl':'ensg',
'Ensembl_PRO':'ensp',
#'GeneID':'ncbig',
#'Gene_Name':'gnm',
#'Gene_Synonym':'gsnm',
#'NCBI_TaxID':'txn',
'RefSeq':'rfsq', # only proteins
'UniParc':'uparc',
#'UniProtKB-ID':'upid',
}
tsvkeys = {
0: {'upac': ';'},
1: {'upid': ';'},
2: {'gnm': ';'},
3: {'gsnm': ';'},
4: {'spnm': ';'},
5: {'txn': ';'},
6: {'def': ''},
7: {'pome': ';'},
8: {'score': ';'},
9: {'pubmed': ';'},
#10: {'ncbig': ';'},
##11: {'rfsq': ';'}, # e.g. NP_001257413.1 [O15178-2];
##12: {'enst': ';'}, # e.g. ENST00000296946 [O15178-1];
}
dirpth = sys.argv[1]
tsv = sys.argv[2]
idm = sys.argv[3]
idms = glob.glob(dirpth + '*.' + idm)
tsvs = glob.glob(dirpth + '*.' + tsv)
mapA = {}
print(time.ctime(time.time()))
if len(sys.argv) > 4:
    mapApth = sys.argv[4]
    mapA = Z.xmap(mapApth, 0, 2) # RefSeq AC=>GN
    count = len(mapA.keys())
    print(str(mapApth) + ':' + str(count))

for idmpth in idms:
    idmdat = Z.idmap(idmpth, 0, 2, idmkeys, mapA)
    count = len(idmdat.keys())
    if not count:
        print('WARNING: prot.py: no data in:' + idmpth)
        continue
    print(str(idmpth) + ':' + str(count))
    ## pprint.pprint(idmdat['ENSG00000156113']) # 70 accs

    tsvpth = idmpth.replace(idm, tsv)
    tsvdat = Z.xflex(tsvpth, '\t', tsvkeys, 0) # map not needed
    count = len(tsvdat.keys())
    if not count: # never 0 ?? TODO
        print('WARNING: prot.py: no data in:' + tsvpth)
    print(str(tsvpth) + ':' + str(count))

    outpth = idmpth.replace(idm, 'nt')
    bgwmap = dat2rdf(idmdat, tsvdat, outpth ) # Attn: tsvdat going to be modified !!! ???
    for key in bgwmap.keys():
        count = len(bgwmap[key].keys())
        print('bgwmap:' + key + ':' + str(count))
    xrfpth = idmpth.replace(idm, 'xrf')
    with open(xrfpth, 'w') as fp:
        json.dump(bgwmap, fp, indent=2)
    fp.close()

print(time.ctime(time.time()))
