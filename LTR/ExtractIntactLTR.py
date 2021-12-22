#!/usr/bin/env python3
#Extract intact LTR sequence based on GFF file
#args[1]: pass.list.gff3 args[2]: genome.fasta args[3]:output.fasta
#AUTHOR:GAOCY
import sys

args=sys.argv


gff_list=[]
with open(args[1],'r') as gff:
    for line in gff:
        lin=line.strip().split('\t')
        if len(lin)!=1 and lin[2]=='repeat_region':
            gff_list.append((lin[0],lin[3],lin[4]))
            
fa_dic={}
with open(args[2],'r') as fa:
    for s in fa:
        s=s.strip()
        if s.startswith('>'):
            sname=s[1:]
            fa_dic[sname]=''
        else:
            fa_dic[sname]+=s
with open(args[3],'w') as w:
    for info in gff_list:
        chr_=info[0]
        start=info[1]
        end=info[2]
        name=f'{chr_}_{start}-{end}'
        fasta=fa_dic[chr_][int(start)-1:int(end)]
        w.write('>'+name+'\n')
        w.write(fasta+'\n')
