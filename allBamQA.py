from Bio import SeqIO
import sys
sys.path.remove('/pub/tool/RSeQC-2.6.1/lib64/python2.7/site-packages/RSeQC-2.6.1-py2.7-linux-x86_64.egg')
import pysam

samfile = pysam.AlignmentFile(sys.argv[3], "r" )

def rateMap(l):
 d={}
 for i in l:
  k,v=i[0],i[1]#499S
  if k not in d.keys():
   d[k]=v
  else:
   d[k]+=v
 if 1 in d.keys():#1:insert 2:deletion 4:soft clip 0:M
  readlength=d[0]+d[1]+d[4]#M+I+S=1313
  maprate="%0.4f"%((d[1]+d[0])/readlength)#I+M
  if 2 in d.keys():
   indel_rate="%0.4f"%((d[1]+d[2])/d[0])
  else:
   indel_rate="%0.4f"%(d[1]/d[0])
 else:
  readlength=d[0]+d[4]#M+S
  maprate="%0.4f"%(d[0]/readlength)#M
  if 2 in d.keys():
   indel_rate="%0.4f"%(d[2]/d[0])
  else:
   indel_rate=0
 return maprate,readlength,indel_rate,d[0]#d[0]==M

refs={}
#fastq=
def fastAQ2dict(f,syb):
 if syb=="fastq":
  seq_records = SeqIO.parse(f,"fastq")
 else:
  seq_records = SeqIO.parse(f,"fasta") 
 d={}
 #seq_records_dict={d[seq_record.id]=seq_record.seq for seq_record in seq_records}
 for seq in seq_records:
  d[seq.id]=seq.seq
 return d

dref=fastAQ2dict(sys.argv[1],"fasta")
dfastq=fastAQ2dict(sys.argv[2],"fastq")

def matchRefQuery(q,r,l):
 match=0;mismatch=0
 for poses in l:
  qpos,rpos=poses[0],poses[1]
  qgene,rgene=q[qpos-1],r[rpos-1]
  if qgene==rgene:
   match+=1
  else:
   mismatch+=1 
 # print(q[qpos-1],r[rpos-1]) 
 return(match,mismatch)

for read in samfile:
 map_rate,readlength,indel_rate,M=rateMap(read.cigar)

 query_seq=dfastq[read.query_name]
 ref_seq=dref[read.reference_name]
 #qlength=read.query_length
 listpos=read.get_aligned_pairs(matches_only=True) 
 match,mismatch=matchRefQuery(query_seq,ref_seq,listpos)#get the number of match and mismatch
 
 distribute1,distribute2=read.reference_start,read.reference_end
 #map_rate,readlength,indel_rate,M=rateMap(read.cigar)
 match_rate="%0.4f"%(match/M)
 mismatch_rate="%0.4f"%(mismatch/M)
 print(read.query_name, read.query_length,read.reference_name,read.reference_length,distribute1,distribute2,map_rate,match_rate,mismatch_rate,indel_rate)
 #print(read.query_name, read.query_length, read.reference_name, read.reference_length,read.template_length,read.query_alignment_length,read.query_alignment_start,read.query_alignment_end,read.reference_start,read.reference_end,read.cigar)
 #rate,readlength=rateMap(read.cigar)
 #print(read.query_name,read.query_alignment_length,read.reference_name,read.reference_length,rate,readlength)
 #print(read.query_name,read.query_alignment_length,read.reference_name,read.reference_length,len(read.get_reference_positions()))
 ######get reference sequence
 #print(read.query_name,read.reference_name,len(read.query_sequence),len(read.query_alignment_sequence),read.get_reference_sequence())
 ######get the match point
 #print(read.query_name,read.reference_name,len(read.query_sequence),len(read.query_alignment_sequence),len(read.get_aligned_pairs(matches_only=True)),read.get_aligned_pairs(matches_only=True))


