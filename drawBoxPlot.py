import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns

def newdf(fil):
 df=pd.read_table(fil,sep=" ",names=['query_name', 'query_length','reference_name','reference_length','dis1','dis2','map_rate','match_rate','mismatch_rate','indel_rate'])
 diffvalue=(df['query_length'].max()-df['query_length'].min())//5
 minlength=df['query_length'].min()
 maxlength=df['query_length'].max()
 v1=minlength+diffvalue;syb1=str(minlength)+"-"+str(v1)
 v2=minlength+diffvalue*2;syb2=str(v1+1)+"-"+str(v2)
 v3=minlength+diffvalue*3;syb3=str(v2+1)+"-"+str(v3)
 v4=minlength+diffvalue*4;syb4=str(v3+1)+"-"+str(v4)
 syb5=str(v4+1)+"-"+str(maxlength)
 def distinguishSYB(length):
  if length >= minlength and length <= v1:
   return syb1
  elif length > v1 and length <= v2:
   return syb2
  elif length > v2 and length <= v3:
   return syb3
  elif length > v3 and length <= v4:
   return syb4
  elif length > v4 and length <= maxlength:
   return syb5
 df['lengthSyb']=df['query_length'].apply(distinguishSYB)
 p1=df['lengthSyb'][df['lengthSyb']==syb1].count()
 p2=df['lengthSyb'][df['lengthSyb']==syb2].count()
 p3=df['lengthSyb'][df['lengthSyb']==syb3].count()
 p4=df['lengthSyb'][df['lengthSyb']==syb4].count()
 p5=df['lengthSyb'][df['lengthSyb']==syb5].count()
 return df,p1,p2,p3,p4,p5

df,p1,p2,p3,p4,p5=newdf(sys.argv[1])
name=(sys.argv[1]).split(".")[0]
print(name,p1,p2,p3,p4,p5)
plotrate="indel_rate"
ax=sns.boxplot(x="lengthSyb",y=plotrate,hue="reference_name",data=df,palette="Set3")
plt.savefig(name+"_"+plotrate+".pdf",dpi=200)
