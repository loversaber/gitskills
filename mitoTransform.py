from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

rcrs=SeqIO.read("/home/liuqi/work/reference/mito_38.fa","fasta")

#mutrcrs=MutableSeq(str(rcrs.seq))
#doublemito=Seq(str(mutrcrs)+str(mutseq))
doublemito=rcrs.seq+rcrs.seq

def mitoPos(a,b):
 rec=doublemito[a-500:b+500]
 return rec,a,b

r2,r21,r22=mitoPos(14258-500,16569+14276+500)
r3,r31,r32=mitoPos(10858-500,16569+10653+500)
print(len(r2),len(r3))

recs=(rcrs,SeqRecord(r2,id="MT_r2",description=str(r21+500)+"f_"+str(r22-16569-500)+"r"+"_"+str(len(r2))),SeqRecord(r3,id="MT_r3",description=str(r31+500)+"f_"+str(r32-16569-500)+"r"+"_"+str(len(r3))))

otfil=open("threeMito.fa","w")
#rec2=SeqRecord(r2,id="MT_r2",description=str(r21+500)+"f_"+str(r22-16569-500)+"r")
SeqIO.write(recs,otfil,"fasta")

