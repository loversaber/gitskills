import pandas as pd
import os
from sys import argv
from functools import reduce

def depthMTMean(f):#/home/liuqi/work/4generations/result/test/mem_p_2/45A_DA_BS_2.depth
 print("depthMTMean",f)
 name_list=f.split(".")
 print(name_list)
 name="_".join((name_list[0].split("/")[-1]).split("_")[0:3])#'45C_GM_BL'
 print(name)
 df=pd.read_table(f,"\t",header=None,names=["ref","loc","depth"],index_col=1)
 if "mem" not in name_list[1:]:####should juge bam file
  return name,df['depth'].mean(axis=0)
 else:
  print(len(df.index),"0")
  l1=list(range(1,501));l2=list(range(16570,17070))
  match_pos=dict(zip(l1,l2)) 
  for key in match_pos.keys():
   df.loc[key,'depth']+=df.loc[match_pos[key],'depth']
   df.drop(labels=match_pos[key],axis=0,inplace=True)
  print(len(df.index),"1")
  return name,df['depth'].mean(axis=0)

def depthGet(f):#fil.endswith(".bam")
 print(f)
 name=f.split(".")[0]
 print(name)
 name2=f.split(".")[1:]
 #if re.search("sorted",name2)!=None:
 if "sorted" not in name2:
  cmd0="samtools sort -@ 40 -o %s %s"
  os.system(cmd0%(f+".sorted.bam",f))
  cmd1="samtools index -@ 40 %s"
  os.system(cmd1%(f+".sorted.bam"))
  cmd2="samtools depth -a -r MT %s > %s"
  os.system(cmd2%(f+".sorted.bam",name+".depth"))
  os.remove(f+".sorted.bam")
  os.remove(f+".sorted.bam.bai")
 else:
  cmd1="samtools index %s"
  os.system(cmd1%(f))
  cmd2="samtools depth -a -r MT %s > %s"
  os.system(cmd2%(f,name+".depth"))
  os.remove(f+".bai")

 depthfil=name+".depth"
 df=pd.read_table(depthfil,"\t",header=None,names=["ref","loc","depth"],index_col=1)
 nameend="_".join((name.split("/")[-1]).split("_")[0:3])
 if "mem" not in name2:####should juge bam file
  return nameend,df['depth'].mean(axis=0)
 else:
  #print(len(df.index),"0")
  l1=list(range(1,501));l2=list(range(16570,17070))
  match_pos=dict(zip(l1,l2))
  for key in match_pos.keys():
   df.loc[key,'depth']+=df.loc[match_pos[key],'depth']
   df.drop(labels=match_pos[key],axis=0,inplace=True)
  #print(len(df.index),"1")
  return nameend,df['depth'].mean(axis=0) 

def pathDF(path):
 d={}
 lfils=[fil for fil in os.listdir(path) if fil.endswith(".bam")]
 depth_syb=path.split("/")[-1]
 #depth_syb="_".join(lfils[0].split(".")[1:])
 #name=fil.split(".")[0]
 for fil in lfils:
  name,depthmean=depthGet(path+'/'+fil)
  d[name]={}
  d[name][depth_syb]=depthmean
 df=pd.DataFrame(d).T
 return df

def allHeter():
 path="/home/liuqi/work/4generations/result/filtered_xls"
 xlsfils=[i for i in os.listdir(path) if i.endswith("_2_plus_filtered.xls")]
 d={}
 for fil in xlsfils:
  name="_".join(fil.split("_")[0:3])
  df=pd.read_table(path+"/"+fil,sep="\t",header=0)
  nheter=(df['Hmarker']==1).count()
  d[name]={}
  d[name]['nheter']=nheter
 #d['22D_DA_BS']['nheter']=0
 #d['40D_DA_BS']['nheter']=0# not work
 df=pd.DataFrame(d).T
 df.loc['22D_DA_BS']=0
 df.loc['40D_DA_BS']=0
 return(df)

with open(argv[1],"r")as paths:
 dl=[]
 for path in paths: 
  path=path.strip()
  dl.append(pathDF(path))
  #print(pathDF(path))
 dfheter=allHeter()
 dl.append(dfheter)
 #print(dl)
 df_final = reduce(lambda left,right: pd.merge(left,right,left_index=True,right_index=True), dl)
 #print(df_final)
 df_final = df_final.sort_index(axis=1)
 df_final.to_csv("depthallfile.txt",sep="\t",index=True,float_format="%.4f")
