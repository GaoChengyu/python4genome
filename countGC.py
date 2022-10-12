from pyfasta import Fasta
import argparse

def getFasta(fafile):
    fa=Fasta(fafile)
    return fa

def getBin(fafile,binLen:int,stepLen:int):
    fa=getFasta(fafile)
    
    for seqid in fa.keys():
        start=0
        falen=len(fa[seqid])
        end=start+binLen
        while start < falen:
            binSeq=fa[seqid][start:end]
            if end >falen:
                end= falen
            range_len=f"{start+1}-{end}"
            start+=stepLen
            end=start+binLen
            yield seqid,binSeq,range_len

def getGCinBin(fafile,binLen:int,stepLen:int):
    for seqid,binSeq,range_len in getBin(fafile,binLen,stepLen):
        seqcount=len(binSeq)
        GCcount=(binSeq.count("G")+binSeq.count("C"))/seqcount
        yield seqid,GCcount,range_len


def getGCinGenome(fafile):
    fa=getFasta(fafile)
    allseqlen=0
    allGClen=0
    for key in fa.keys():
        seqlen=len(fa[key])
        GCnum=str(fa[key]).count("G")+str(fa[key]).count("C")
        allseqlen+=seqlen
        allGClen+=GCnum
    allGCcount=allGClen/allseqlen
    return allGCcount

def getGCperChr(fafile):
    fa=getFasta(fafile)
    for key in fa.keys():
        seqlen=len(fa[key])
        GCcount=(str(fa[key]).count("G")+str(fa[key]).count("C"))/seqlen
        yield key,GCcount

def output(fafile,binLen:int,stepLen:int,path:str):
    with open(f'{path}/GCperBin.txt','w') as w1:
        for seqid,GCcount,range_len in getGCinBin(fafile,binLen,stepLen):
            w1.write(f'{seqid}\t{GCcount}\t{range_len}\n')
    with open(f'{path}/GC_of_all.txt','w') as w2:
        GCcount=getGCinGenome(fafile)
        w2.write(f'{GCcount}')
    with open(f'{path}/GCperChr.txt','w') as w3:
        for seqid,GCcount in getGCperChr(fafile):
            w3.write(f'{seqid}\t{GCcount}\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('inputFile',type=str,help='输入fasta文件')
    parser.add_argument('outputFilePath',type=str,help='输出文件所在目录')
    parser.add_argument('--bin',type=int,default='500',help='选择滑窗范围')
    parser.add_argument('--step',type=int,default='500',help='选择步长，即每次滑窗移动的距离')
    args=parser.parse_args()
    output(args.inputFile,args.bin,args.step,args.outputFilePath)


if __name__=='__main__':
    main()

