#!/usr/bin/env python
#coding: utf8


####import module
from __future__ import division
import sys, os
import re
import glob
import gzip
import json
import codecs
import string
import random
import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import seaborn as sns
#####Document Description
''' This script is created for calculating block Q30, and mapping rate
'''
##########log#############

Date = '2018-06-05'
Version='0.3.0'

usage = '''
      Version %s Ya Ding %s

      Usage: %s <samfile.sam> <version> <Output Prefix>>STDOUT
''' %(Version, Date, os.path.basename(sys.argv[0]))

##########################change log#################
#v0.2 able to calulate BlockQ30,mappingRate and cycle Q30, two report generated
#v0.3 final version: add GC content and errorRate
#rename v0.3

def block_find(ID):
    a=[]
    b=[]
    blk=0
    for i in range(0,8):
        for j in range(0,8):
            a.append(y[i]*x[j])

    for i in range(1,65):
            b.append(sum(a[:i])-1)

    if ID <= x[0]*y[0]-1:
        blk= 1
    else:
        for q in range(0,63):
            if ID<=b[q+1] and ID>b[q]:
                blk= q+2
    return blk


def Q30(q_str):
    q=0
    for i in range(len(q_str)):
        if ord(q_str[i])-33>30:
            q += 1
    return q

def Q30_base(base):
    q=0
    if ord(base)-33>30:
        q+=1
    return q

def mapIdx(num):
    unmaped = [73,133,181,89,117,153,77,69,137,141]
    if num <=16:
        if num !=0 and num !=16:
            return 0
        else:
            return 1
    else:
        if num not in unmaped:
            return 1
        else:
            return 0

def GC_cal(p_str):
    gc=0
    for i in range(len(p_str)):
        if p_str[i] == 'G' or p_str[i] == 'C':
            gc += 1
        else:
            continue
    return gc

def Err_cal(e_str):
    err = 0
    err_str = e_str.split(':')[-1]
    for i in range(len(err_str)):
        if err_str[i].isalpha():
            err +=1
        else:
            continue
    return err

def bam_read(bamfile):
    FOV={} # store the number of reads
    Qval={}
    bl=[[] for i in range(65)]
    html_tmp={}         # for hmtl template
    mapR={}
    cycleQ={}
    GC={}
    Err={}

    with open (bamfile, 'r') as bam:
        bam.readline()
        bam.readline()
        while True:
            idline= bam.readline()
            if len(idline.split("\t")) <10:
                break
            idfor = idline.split("\t")[0]
            q_str=idline.split("\t")[10]
            p_str=idline.split("\t")[9]
            flag=int(idline.split("\t")[1])
            if mapIdx(flag):
                for i in range(len(idline.split("\t"))):
                    if 'MD:Z' in idline.split("\t")[i]:
                        e_str = idline.split("\t")[i]
                    else:
                        continue
            length=len(idline.split("\t")[9])
            m = re.search('C\d{3}R\d{3}', idfor)
            idx = m.start()
            fov = idfor[idx:idx+8]
            if '/' in idfor:
                rid = int(idfor.split('/')[0][idx+8:])
            elif idfor[idx+8]=='_':
                rid = int(idfor.split('_')[-1])
            else:
                rid = int(idfor.strip()[idx+8:])
            if rid <Total-1:
                if fov not in FOV.keys():
                    FOV[fov]={}
                    Qval[fov]={}
                    mapR[fov]={}
                    cycleQ[fov]={}
                    GC[fov]={}
                    Err[fov] ={}
                    for j in range(length):
                        cycleQ[fov][j+1]={}
                        for i in range (1,65):
                            cycleQ[fov][j+1][i]=0


                    for i in range (1,65):
                        FOV[fov][i]=0
                        Qval[fov][i]=0
                        mapR[fov][i]=0
                        GC[fov][i]=0
                        Err[fov][i] = 0
                    blk= block_find(rid)
                    mapidx =mapIdx(flag)

                    FOV[fov][blk] = FOV[fov][blk] +1   # count total reads
                    mapR[fov][blk] = mapR[fov][blk] + mapidx  # count mapped reads

                    if mapidx == 1:
                        error = Err_cal(e_str)
                        gc = GC_cal(p_str)
                        qval = Q30(q_str)
                        Qval[fov][blk] = Qval[fov][blk] + qval
                        GC[fov][blk] = GC[fov][blk] +gc
                        Err[fov][blk] = Err[fov][blk] + error
                        for i in range(length):
                            p = Q30_base(q_str[i])
                            cycleQ[fov][i+1][blk]=cycleQ[fov][i+1][blk]+p
                    else:
                        continue


                else:
                    blk= block_find(rid)
                    mapidx =mapIdx(flag)
                    FOV[fov][blk] = FOV[fov][blk] +1
                    mapR[fov][blk] = mapR[fov][blk] + mapidx
                    if mapidx == 1:
                        qval = Q30(q_str)
                        gc = GC_cal(p_str)
                        error = Err_cal(e_str)
                        FOV[fov][blk] = FOV[fov][blk] +1
                        Qval[fov][blk] = Qval[fov][blk] + qval
                        mapR[fov][blk] = mapR[fov][blk] + mapidx
                        GC[fov][blk] = GC[fov][blk] +gc
                        Err[fov][blk] = Err[fov][blk] + error
                        for i in range(length):
                            p = Q30_base(q_str[i])
                            cycleQ[fov][i+1][blk]=cycleQ[fov][i+1][blk]+p
                    else:
                        continue
            else:
                continue
    keylist = sorted(list(Qval.keys()))

    for fov in keylist:

        for j in range(length):
            for i in range(1,65):
                if cycleQ[fov][j+1][i] != 0:
                    cycleQ[fov][j+1][i]=cycleQ[fov][j+1][i]/mapR[fov][i]


        for i in range(1,65):
            if Qval[fov][i] != 0:

                Qval[fov][i] = Qval[fov][i] / (mapR[fov][i]*length)

                GC[fov][i] = GC[fov][i]/(mapR[fov][i]*length)

                Err[fov][i] = Err[fov][i]/(mapR[fov][i]*length)

                mapR[fov][i] = mapR[fov][i] / FOV[fov][i]
            else:
                continue


    return Qval,mapR, keylist,cycleQ,length,GC,Err

def html_gen(Qval,mapR,keylist,GC,Err):
    html_tmp={}
    CR_list=[]
    for i in range(1,9):
        for j in range(1,9):
            CR_list.append('C0{}R0{}'.format(i,j))

    for fov in keylist:
        html_tmp[fov]=[[] for i in range(64)]
        for h in range(64):
            html_tmp[fov][h].append(CR_list[h])
            html_tmp[fov][h].append(Qval[fov][h+1])
            html_tmp[fov][h].append(mapR[fov][h+1])
            html_tmp[fov][h].append(GC[fov][h+1])
            html_tmp[fov][h].append(Err[fov][h+1])
            #html_tmp[fov][h].append(random.uniform(0,1))
            #html_tmp[fov][h].append(random.uniform(0,1))

    return html_tmp

def html2_gen(cycleQ,keylist,length):
    html2_tmp={}
    CR_list=[]
    for i in range(1,9):
        for j in range(1,9):
            CR_list.append('C0{}R0{}'.format(i,j))

    for fov in keylist:
        html2_tmp[fov]=[[] for i in range(64)]
        for h in range(64):
            html2_tmp[fov][h].append(CR_list[h])
            for j in range(length):
                html2_tmp[fov][h].append(cycleQ[fov][j+1][h+1])

    return html2_tmp

def main():
    if len(sys.argv) !=4:
        print (usage)
        sys.exit(1)

    bamfile=sys.argv[1]
    outfile=sys.argv[3]
    bamdir = os.path.dirname(bamfile)
    pyfile = os.path.dirname(sys.argv[0]) 
    
    global x, y,Total
    if sys.argv[2] == 'V2':
        x=[67,109,165,193,193,165,109,67]
        y=[45,61,125,173,173,125,61,45]
        Total =862943
    elif sys.argv[2] == 'V01':
        x=[73,93,107,117,117,107,93,73]
        y=[73,93,107,117,117,107,93,73]
        Total =608399
    elif sys.argv[2] == 'V1':
        x=[62,79,92,102,102,92,79,62]
        y=[52,60,74,87,87,74,60,52]
        Total =365819
    elif sys.argv[2] == 'V02':
        x=[57,102,132,192,192,132,102,57]
        y=[45,81,105,153,153,105,81,45]
        Total =741887
    elif sys.argv[2] == 'V8':
        x=[73,130,149,225,225,149,130,73]
        y=[73,130,149,225,225,149,130,73]
        Total =1331715
    else:
        x=[189,237,285,309,309,285,237,189]
        y=[189,237,285,309,309,285,237,189]
        Total =4161599

    os.chdir(bamdir)
    Qval, mapR, fov_list, cycleQ, length,GC,Err = bam_read(bamfile)
    cycle_block=html2_gen(cycleQ,fov_list,length)

    Fov_block=html_gen(Qval,mapR,fov_list,GC,Err)
    option1 = ['BlockQ30','MappingRate','GC content','ErroRate']
    jsonData1 = {}
    jsonData1['blockData'] = json.dumps(Fov_block)
    jsonData1['option'] = json.dumps(option1)
    fhIn1 =codecs.open('{}\heatmap_template.html'.format(pyfile),encoding='utf-8')
    temp1 = fhIn1.read()
    newStr1 = string.Template(temp1)
    newStr1 = newStr1.safe_substitute(jsonData1)
    if os.path.exists(outfile):
        os.remove(outfile)
    fhout = codecs.open('{}.html'.format(outfile), 'w',encoding='utf-8')
    fhout.write(newStr1)
    fhout.close()


    ####second template for cycle_blockQ30################
    option=[]
    for i in range(length):
        option.append('S'+str(i+1))
    jsonData2={}
    jsonData2['blockData'] = json.dumps(cycle_block)
    jsonData2['option'] = json.dumps(option)
    fhIn2 =codecs.open('{}\heatmap_template2.html'.format(pyfile),encoding='utf-8')
    temp2 = fhIn2.read()
    newStr2 = string.Template(temp2)
    newStr2 = newStr2.safe_substitute(jsonData2)
    outfile2='{}_cycleQ.html'.format(outfile)
    if os.path.exists(outfile):
        os.remove(outfile)
    fhout = codecs.open(outfile2, 'w',encoding='utf-8')
    fhout.write(newStr2)
    fhout.close()






if __name__ == '__main__':
    main()
