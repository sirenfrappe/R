#!/bin/bash -e

echo `date`

fasta=$1
id=${fasta##*/}
id=${id%.*}
path_qianzhui=$(pwd)

#SPIDER3
echo SPIDER3 started
/mnt/d/linux/SPIDER3-numpy-server/run_list.sh $fasta
cp /mnt/d/linux/$id.fasta.spd33 /mnt/d/linux/DeepConPred2/data/$id.spd33
rm $id*
echo SPIDER3 finished

#CCMPred
#generate alignments
echo CCMPred started
#$path_qianzhui/scripts/generate-alignments.pl $fasta $path_qianzhui/scripts/out/$id
#/mnt/d/DNCON2/CCMpred/bin/ccmpred -t 8 $path_qianzhui/scripts/out/$id/$id.aln $path_qianzhui/scripts/out/$id/$id.ccmpred > $path_qianzhui/scripts/out/$id/$id.log
#cp $path_qianzhui/scripts/out/$id/$id.ccmpred $path_qianzhui/DeepConPred2/data
echo CCMPred finished

#psiblast
echo psiblast started
#/public/home/SZUYanwy/DNCON2/ncbi-blast-2.2.25+/bin/psiblast -query $fasta -evalue .001 -inclusion_ethresh .002 -db $path_qianzhui/../DNCON2/databases/nr90-2012/nr90 -num_iterations 3 -num_threads 8 -seg yes -out_ascii_pssm $path_qianzhui/DeepConPred2/data/$id.PSSM
echo psiblast finished

#deepConPred
echo DeepConPred2 started
cd $path_qianzhui/DeepConPred2
#fasta
cp ./data/$id.fasta ./data/$id 
sed  -i '1d' ./data/$id 
sed -i ':a;N;s/\n//g;ta' ./data/$id
sed -i '1i >'"$id"' ' ./data/$id
#spd33
cp ./data/$id.spd33 ./data/$id.spd33_old
Rscript /mnt/d/linux/scripts/transSpd33.r ./data/$id.spd33
sed -i '1i # SEQ SS ASA Phi Psi Theta(i-1=>i+1) Tau(i-2=>i+2) HSE_alpha_up HSE_alpha_down P(C) P(H) P(E)' ./data/$id.spd33
python $path_qianzhui/DeepConPred2/DeepConPred2.py $id
echo DeepConPred2 finished


echo $id finished
echo `date`
