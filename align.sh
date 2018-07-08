minimap2=/home/liuqi/tool/minimap2-2.11_x64-linux/minimap2

blasr=/home/liuqi/miniconda3/bin/blasr

fastq=/home/liuqi/tmp/nanopore
#ref=/home/liuqi/work/reference
parafly=/home/liuqi/tool/parafly-r2013-01-21/bin/bin/ParaFly
ref3=/home/liuqi/tmp/nanopore

#$minimap2 -d $ref/mito_38.mmi $ref/mito_38.fa

#$minimap2 -d $ref3/threeMito.mmi $ref3/threeMito.fa
for i in `find $fastq -size +60M -name *.fastq`
do
a=${i%%.fastq}
name=`echo $a|awk -F"/" '{print $NF"_"$(NF-1)}'`
b=${a##*/}
echo $name
###minimap2 
#echo "$minimap2 -ax map-ont $ref3/threeMito.fa $i > $name".sam"" >> map.txt

###blasr 
echo "$blasr $i $ref3/threeMito.fa --bam --out $b"_blasr.bam" --nproc 40" >> map_blasr.txt

#ls -alth $i
done

#$parafly -c map.txt -CPU 10
$parafly -c map_blasr.txt -CPU 10
