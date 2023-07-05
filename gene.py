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
    'gn',
    ]
    ### Object properties ###
    props = [
    'gn2gp',
    'gn2txn',
    'sth2ori', # TODO
    'sth2evd',
    'sth2src',
    'sth2clm',
    'mbr2lst',
    'cls2prn'
    ]
    ### Annotation properties ###
    aprops = [
    'sth2nm',
    'sth2syn',
    'sth2dfn',
    ]
    itemUs = {'self': '<http://rdf.biogateway.eu/graph/gene>'}
    buff = Z.header(itemUs, props, aprops, olders)
    fp = open(outpth, 'w')
    fp.write(buff)
    buff = ''
    ###########################################################################
    ### cardinality accession:GeneNmae is N:M !!!

    prnU = itemUs['gn']
    bagU = '<' + Z.uris['rdf'] + 'Bag' + '>'
    txnBU = Z.uris['ncbitx']
    geneBU = Z.uris['gene']
    protBU = Z.uris['prot']
    #prBU = Z.uris['pr'] # do NOT delete!
    mysrc = 'uniprot'
    srcBU = Z.uris[mysrc] # do NOT delete!
    srcU = '<' + srcBU + '>'
    outcnt = 0
    # datafields to be treated as xrfs
    insts = (
    'ensg',
    'ncbig',
    )
    # data fields to be transferred to isoforms
    onenmsp = 'gnm'
    bgwmap = {'bgw': {}}

    buff += itemUs['self'] + ' ' + itemUs['sth2src'] + ' ' + srcU + Z.END
    idminv = Z.upac2bac(idmdat)
    for onegnm in sorted(tsvdat.keys()):
        unidat = tsvdat[onegnm]
        onebacs = sorted(unidat['upac'].keys())
        for onebac in onebacs:
            # propagating data to isoforms
            try:
                idmbac = idminv[onebac]
            except KeyError:
                continue
            upacs = sorted(idmbac.keys())
            for upac in upacs:
                if upac not in unidat['upac'].keys():
                    unidat['upac'][upac] = 1
                idmac = idmbac[upac]
                for mynmsp in idmac.keys():
                    if mynmsp not in unidat.keys():
                        unidat[mynmsp] = {}
                    myids = idmac[mynmsp].keys()
                    for myid in myids:
                        unidat[mynmsp][myid] = 1

        ### sanity checks and filtering
        txns = sorted(unidat['txn'].keys())
        if len(txns) != 1: # couunt:0
            Z.mapAlert(txns, 'txns', onegnm, onenmsp)
            continue
        mynmsp = 'uparc'
        try:
            gpids = sorted(unidat[mynmsp].keys())
        except KeyError: # 5985 gene names
            #print('data2rdf:KeyError:unidat[' + mynmsp + ']:' + onegnm)
            continue # effecively filtering by RefProt
        mynmsp = 'spnm'
        if mynmsp in unidat.keys():
            spnms = sorted(unidat[mynmsp].keys())
        if len(spnms) == 0: # sic - for some no
            continue
        if len(spnms) > 1: # none
            Z.mapAlert(spnms, 'spnms', onegnm, onenmsp)
        mynmsp = 'pome'
        if mynmsp in unidat.keys():
            pomes = sorted(unidat[mynmsp].keys())
        if len(pomes) != 1:
            # none: 2, mutiple: 97; all unique
            Z.mapAlert(pomes, 'pomes', onegnm, onenmsp)
        somes = {}
        for pome in pomes:
            bits = pome.split(' ')
            if len(bits) != 3:
                continue # skipping 'Unplaced' ?? TODO
            mypome = bits[0]
            some = bits[2]
            if mypome not in somes.keys():
                somes[mypome] = {}
            somes[mypome][some] = 1
        allpomes = sorted(somes.keys())
        if len(allpomes) != 1: # 790, all none and unique
            #Z.mapAlert(allpomes, 'allpomes', onegnm, onenmsp)
            continue # seem 'Unplaced' TODO
        onesome = somes[allpomes[0]]

        ### exporting
        ## defining properties
        txn = txns[0]
        for mysome in onesome.keys():
            if len(mysome) == 1:
                smnm = 'chr-0' + mysome
            else:
                smnm = 'chr-' + mysome
            oneid = txn + '/' + smnm + '/' + onegnm + '/' # terminating with '/' for bags
            oneU = '<' + geneBU + oneid + '>'
            mypdcU = itemUs['cls2prn']
            #buff += oneU + ' ' + mypdcU + ' ' + prnU + Z.END
            #outcnt += 1
            buff += oneU + ' ' + mypdcU + ' ' + bagU + Z.END
            outcnt += 1
            # taxon
            mynmsp = 'txn'
            mypdcU = itemUs['gn2txn']
            myBU = txnBU
            out = Z.me2you(unidat, mynmsp, oneU, mypdcU, myBU)
            buff += out[0]
            outcnt += out[1]
            # translation products
            mynmsp = 'uparc'
            mypdcU = itemUs['gn2gp']
            for gpid in  gpids:
                gpid = oneid + gpid
                gpU = '<' + Z.uris['prot'] + gpid + '>'
                buff += oneU + ' ' + mypdcU + ' ' + gpU + Z.END
                outcnt += 1
    
            ## source
            # me2you() cannot be used here
            for upac in upacs:
                upbac = upac.split('-')[0]
                if upbac != upac:
                    continue # using only canonical ACs
                onesrcU = '<' + srcBU + upbac + '>'
                mypdcU = itemUs['sth2ori']
                buff += oneU + ' ' + mypdcU + ' ' + onesrcU + Z.END
                outcnt += 1
    
            ## name
            mypdcU = itemUs['sth2nm']
            onenm = onegnm
            buff += oneU + ' ' + mypdcU + ' ' + Z.DQ + onenm + Z.DQ + Z.END
            outcnt += 1
    
            ## synonyms
            mypdcU = itemUs['sth2syn']
            mynmsp = 'gsnm'
            out = Z.me2lbl(unidat, mynmsp, oneU, mypdcU)
            if out:
                buff += out[0]
                outcnt += out[1]
    
            ## description
            allsomes = Z.DS.join(somes[allpomes[0]].keys()) # 2 DSs in the output
            spnm = spnms[0]
            dfn = 'Gene ' + onegnm + ' from ' + spnm + ' chromeome '  + allsomes
            mypdcU = itemUs['sth2dfn']
            buff += oneU + ' ' + mypdcU + ' ' + Z.DQ + dfn + Z.DQ + Z.END
            outcnt += 1
    
            ## publications
            mypdcU = itemUs['sth2evd']
            mynmsp = 'pubmed'
            myBU = Z.uris[mynmsp]
            if mynmsp in unidat.keys():
                out = Z.me2you(unidat, mynmsp, oneU, mypdcU, myBU)
            if out:
                buff += out[0]
                outcnt += out[1]
    
            ## xrfs
            mypdcU = itemUs['sth2clm']
            for mynmsp in insts:
                if mynmsp not in unidat.keys():
                    continue
                myBU = Z.uris[mynmsp]
                out = Z.me2you(unidat, mynmsp, oneU, mypdcU, myBU)
                if not out:
                    continue
                buff += out[0]
                outcnt += out[1]
                # building bgwmap
                if oneid not in bgwmap['bgw'].keys():
                    bgwmap['bgw'][oneid] = {}
                if mynmsp not in bgwmap.keys():
                    bgwmap[mynmsp] = {}
                if mynmsp not in bgwmap['bgw'][oneid].keys():
                    bgwmap['bgw'][oneid][mynmsp] = {}
                for myid in unidat[mynmsp]:
                    if myid not in bgwmap[mynmsp].keys():
                        bgwmap[mynmsp][myid] = {}
                        bgwmap[mynmsp][myid]['bgw'] = {}
                    bgwmap[mynmsp][myid]['bgw'][oneid] = 1
                    bgwmap['bgw'][oneid][mynmsp][myid] = 1
    
            # members
            mypdcU = itemUs['mbr2lst']
            for mynmsp in insts:
                if mynmsp not in unidat.keys():
                    continue
                myBU = Z.uris[mynmsp]
                out = Z.you2me(unidat, mynmsp, oneU, mypdcU, myBU)
                if not out:
                    continue
                buff += out[0]
                outcnt += out[1]
    
            fp.write(buff) # per gene
            buff = ''

    fp.close()
    print(outpth + ':' + str(outcnt))

    return(bgwmap)
    ###################### END OF dat2rdf #####################################

### Fields
idmkeys = {
'Ensembl':'ensg',
#'Ensembl_PRO':'ensp',
#'GeneID':'ncbig',
#'Gene_Name':'gnm',
#'Gene_Synonym':'gsnm',
#'NCBI_TaxID':'txn',
#'RefSeq':'rfsq', # only proteins
'UniParc':'uparc',
#'UniProtKB-ID':'upid',
}
tsvkeys = {
0: {'upac': ';'},
#1: {'upid': ';'},
2: {'gnm': ';'},
3: {'gsnm': ';'},
4: {'spnm': ';'},
5: {'txn': ';'},
#6: {'def': ''},
7: {'pome': ';'},
#8: {'score': ';'},
9: {'pubmed': ';'},
10: {'ncbig': ';'},
#11: {'rfsq': ';'}, # nucleotides (transxripts?)
#12: {'enst': ';'},
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
    if not count: # this works as expected
        print('WARNING: gene.py: no data in:' + idmpth)
        continue
    print(str(idmpth) + ':' + str(count))
    ## pprint.pprint(idmdat['ENSG00000156113']) # 70 accs

    tsvpth = idmpth.replace(idm, tsv)
    tsvdat = Z.xflex(tsvpth, '\t', tsvkeys, 2)
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
