minimap2=/home/liuqi/tool/minimap2-2.11_x64-linux/minimap2
blasr=/home/liuqi/miniconda3/bin/blasr
fastq=/home/bioinfo/liuqi/tmp/nanopore
parafly=/pub/tool/trinityrnaseq-2.0.6/trinity-plugins/parafly/bin/ParaFly
ref3=/home/bioinfo/liuqi/tmp/nanopore/threeMito.fa

sam=/home/bioinfo/liuqi/tmp/nanopore/minimap2_result
for i in `find $sam -size +60M -name "*.sam"`
do
a=${i%%.sam}
b=${a##*/}
#echo $b
echo "samtools view -h -F 2308 $i > $b"_filter.sam"" >> filter.txt
done

$parafly -c filter.txt -CPU 10
#python3 test/allBamQA.py refFour.fa fastq_runid_61abc24ddb9dd807a5437db7f6632640de638563.fastq lambda_ref4_filter.sam > allreadsRate.xls
for i in `find $fastq -maxdepth 3 -size +60M -name *.fastq`
do
a=${i%%.fastq}
dirname=`echo $a|awk -F"/" '{print $(NF-1)}'`
name_only=`echo $a|awk -F"/" '{print $NF}'`
name=`echo $a|awk -F"/" '{print $NF"_"$(NF-1)}'`
#echo $dirname
#echo $name_only
#ls -alth $sam/$name".sam"
echo "python3 /home/bioinfo/liuqi/tmp/nanopore/pass/barcode06/test/allBamQA.py $ref3 $fastq/pass/$dirname/$name_only".fastq" $sam/$name"_filter.sam" > $name".xls"" >> rateGet.txt
done

$parafly -c rateGet.txt -CPU 10
