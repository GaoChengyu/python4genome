#!/usr/bin/env python3
#author：GAOCHENGYU
#NWAFU

#This script can rename GFF files

import re
import argparse
class RenameGFF():
    def __init__(self,oldgff:str,newgff:str,prefix:str,model:str):
        self.oldgff=oldgff
        self.newgff=newgff
        self.prefix=prefix
        self.model=model

    #读取旧的gff文件
    #Read old GFF files
    def readgff(self):
        oldgff=self.oldgff
        with open(oldgff,'r') as f:
            for line in f:
                yield line

    #形成一个旧基因号对应全部基因信息的字典，键为基因号，值为gene，mRNA，CDS等全部信息
    #注意：此时需要输入一个参数model，作为一个模型识别旧基因号，如旧基因号为FG1G10001,则输入模板为FG
    def processgff(self):
        gffinfo_dic={}
        model=self.model
        for line in self.readgff():
            line=line.strip().split('\t')
            type_=line[2]
            character=line[-1].split(';')
            for name_ in character:
                if type_=='gene':
                    m=re.search(fr'(?<=ID=).*({model}.+)',name_)
                    if m is not None:
                        gffinfo_dic.setdefault(m.group(1),[]).append(line[:-1])
                if type_!='gene':
                    m=re.search(fr'(?<=Parent=).*({model}.+)',name_)
                    if m is not None:
                        gffinfo_dic.setdefault(m.group(1),[]).append(line[:-1])
        return gffinfo_dic

    #根据基因在基因组染色体/scaffold/contig上的位置给基因排序，start越小越靠前
    def sortID(self):
        model=self.model
        IDdic={}
        for line in self.readgff():
            line=line.strip().split('\t')
            chr_=line[0]
            type_=line[2]
            start_=line[3]
            character=line[-1].split(';')
            for name_ in character:
                if type_=='gene':
                    m=re.search(fr'(?<=ID=).*({model}.+)',name_)
                    if m is not None:
                        IDdic.setdefault(chr_,[]).append((start_,m.group(1)))
        sorted_id_dic={}
        for key in IDdic.keys():
            sortlist=sorted(IDdic[key],key=lambda x:int(x[0]),reverse=False)
            sorted_old_id=[i[1] for i in sortlist]
            sorted_id_dic[key]=sorted_old_id
        return sorted_id_dic

    #形成一个旧ID对应新ID的字典
    def rename(self):
        prefix=self.prefix
        sorted_id_dic=self.sortID()
        new_id_dic={}
        for key,value in sorted_id_dic.items():
            chr_num=re.search(r'\d+',key)
            geneNum=0
            for old_id in value:
                geneNum+=1
                geneNum_str=str(geneNum).zfill(4)
                new_id=prefix+chr_num.group()+'G'+geneNum_str+'0'
                new_id_dic[old_id]=new_id
        return new_id_dic

    def make_new_gff(self):
        new_id_dic=self.rename()
        old_gffinfo=self.processgff()
        new_gffinfo=[]
        for key,value in old_gffinfo.items():
            new_id=new_id_dic[key]
            old_info_T=list(zip(*value))
            exon_num_minus=old_info_T[2].count('exon')
            exon_num_plus=1
            cds_num_minus=old_info_T[2].count('CDS')
            cds_num_plus=1
            for info_ in value:
                if info_[2]=='gene':
                    character=f'ID=gene:{new_id};Name={new_id}'
                    info_.append(character)
                    one_line_newinfo='\t'.join(info_)
                    new_gffinfo.append(one_line_newinfo)
                if info_[2]=='mRNA':
                    character=f'ID=mRNA:{new_id};Name={new_id};Parent=gene:{new_id}'
                    info_.append(character)
                    one_line_newinfo='\t'.join(info_)
                    new_gffinfo.append(one_line_newinfo)
                if info_[2]=='exon':
                    if info_[6]=='+':
                        character=f'ID=exon:{new_id}.{str(exon_num_plus)};Parent=mRNA:{new_id}'
                        info_.append(character)
                        one_line_newinfo='\t'.join(info_)
                        new_gffinfo.append(one_line_newinfo)
                        exon_num_plus+=1
                    if info_[6]=='-':
                        character=f'ID=exon:{new_id}.{str(exon_num_minus)};Parent=mRNA:{new_id}'
                        info_.append(character)
                        one_line_newinfo='\t'.join(info_)
                        new_gffinfo.append(one_line_newinfo)
                        exon_num_minus-=1
                if info_[2]=='CDS':
                    if info_[6]=='+':
                        character=f'ID=CDS:{new_id}.{str(cds_num_plus)};Parent=mRNA:{new_id}'
                        info_.append(character)
                        one_line_newinfo='\t'.join(info_)
                        new_gffinfo.append(one_line_newinfo)
                        cds_num_plus+=1
                    if info_[6]=='-':
                        character=f'ID=CDS:{new_id}.{str(cds_num_minus)};Parent=mRNA:{new_id}'
                        info_.append(character)
                        one_line_newinfo='\t'.join(info_)
                        new_gffinfo.append(one_line_newinfo)
                        cds_num_minus-=1
        return new_gffinfo

    def output(self):
        newgff=self.newgff
        with open(newgff,'w') as f:
            new_gffinfo=self.make_new_gff()
            for line in new_gffinfo:
                f.write(line+'\n')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('inputFilePath',type=str,help='输入文件路径')
    parser.add_argument('outputFilePath',type=str,help='输出文件路径')
    parser.add_argument('--profix',type=str,help='新基因名的前缀')
    parser.add_argument('--model',type=str,help='旧基因名的模式前缀')
    args=parser.parse_args()
    renamegff=RenameGFF(args.inputFilePath,args.outputFilePath,args.profix,args.model)
    renamegff.output()


    
if __name__=='__main__':
    main()
