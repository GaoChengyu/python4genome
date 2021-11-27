import re
import argparse


class EXON_CACULATE():
    def __init__(self,gff:str,output_txt:str):
        self.gff=gff
        self.output_txt=output_txt

#读取gff文件为列表
    def read_gff(self):
        gfffile=self.gff
        gfflist=[]
        with open(gfffile,'r') as f:
            for line in f.readlines():
                line= line.strip().split()
                gfflist.append(line)
        return gfflist

#创建一个以id名为key，exon列表为value的字典
    def make_exon_num_dic(self):
        gff=self.read_gff()
        exon_message_dic={}
        for line in gff:
            type_=line[2]
            character=line[-1].split(';')
            if type_=='exon':
                for i in character:
                    m=re.search(r'Parent=.+',i)
                    if m is not None:
                        id_name=m.group()
                        exon_message_dic.setdefault(id_name,[]).append(type_)
        return exon_message_dic

#创建exon数量的列表
    def make_exon_num_list(self):
        exon_message_dic=self.make_exon_num_dic()
        exon_num_list=[]
        for value in exon_message_dic.values():
            exon_num_list.append(str(len(value)))
        return exon_num_list

#计算拥有最多exon基因的exon数
    def find_max_exon_num(self):
        exon_message_dic=self.make_exon_num_dic()
        exon_num_list=[]
        for value in exon_message_dic.values():
            exon_num=len(value)
            exon_num_list.append(exon_num)
        max_num=max(exon_num_list)
        return max_num

#创建字典，key=exon数目，value=具有这个exon数目的基因数
    def calculate_exon_num(self):
        exon_num_list=self.make_exon_num_list()
        max_num=self.find_max_exon_num()
        cal_exon_num_dic={}
        if max_num<=10:
            for site in range(max_num):
                exon_num=exon_num_list.count(str(site))
                cal_exon_num_dic[max_num]=exon_num
        else:
            for site in range(10):
                exon_num=exon_num_list.count(str(site))
                cal_exon_num_dic[str(site+1)]=exon_num
            gt_sum=0
            for str_ in exon_num_list:
                if int(str_)>=10:
                    gt_sum+=1
            cal_exon_num_dic['gt10']=gt_sum
        output_txt=self.output_txt
        with open(output_txt,'w') as f:
            for key in cal_exon_num_dic.keys():
                f.write(str(key)+'\t'+str(cal_exon_num_dic[key])+'\n')

#计算外显子长度分布
    def make_exon_length_dic(self):
        gff=self.read_gff()
        pre_exon_length_dic={}
        for line in gff:
            type_=line[2]
            character=line[-1].split(';')
            exon_length=int(line[4])-int(line[3])
            if type_=='exon':
                for i in character:
                    m=re.search(r'Parent=.+',i)
                    if m is not None:
                        id_name=m.group()
                        pre_exon_length_dic.setdefault(id_name,[]).append(str(exon_length))
        exon_length_dic={}
        for id in pre_exon_length_dic.keys():
            sum_length=0
            for length in pre_exon_length_dic[id]:
                sum_length+=int(length)
            exon_length_dic[id]=str(sum_length)
        return exon_length_dic

    def calculate_exonLength_distribution(self,bin):
        bin=int(bin)
        if 4000//bin==4000/bin:
            block_num=4000//bin
        else:
            block_num=4000//bin+1
        exon_lenth_dic=self.make_exon_length_dic()
        exon_lenth_list=exon_lenth_dic.values()
        exonlengthDistribution_dic={}
        for i in range(block_num):
            start=i*bin+1
            end=(i+1)*bin
            if end>4000:
                end=4000
            range_=str(start)+'-'+str(end)
            length_sum=0
            for length in exon_lenth_list:
                if int(length) >= start and int(length) <= end:
                    length_sum+=1
            exonlengthDistribution_dic[range_]=str(length_sum)
        output_txt=self.output_txt
        with open(output_txt,'w') as f:
            for key in exonlengthDistribution_dic.keys():
                f.write(str(key)+'\t'+str(exonlengthDistribution_dic[key])+'\n')



        


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('inputFilePath',type=str,help='输入文件路径')
    parser.add_argument('outputFilePath',type=str,help='输出文件路径')
    parser.add_argument('--model',type=str,help='选择计算模式 --model CEN 计算外显子数量；--model CED 计算外显子分布')
    parser.add_argument('--bin',type=str,default='500',help='选择滑窗范围')
    args=parser.parse_args()

    if args.model=='CEN':
        cal_exon_num=EXON_CACULATE(args.inputFilePath,args.outputFilePath)
        cal_exon_num.calculate_exon_num()
    if args.model=='CED':
        cal_exonLength_distri=EXON_CACULATE(args.inputFilePath,args.outputFilePath)
        cal_exonLength_distri.calculate_exonLength_distribution(args.bin)


if __name__=='__main__':
    main()
