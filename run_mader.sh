

designs="ER HER2 TNBC relapse relapse_ER relapse_HER2 relapse_TNBC"

#epcy kal2mat -d ./data/design/mader/all/design.tsv -o ./data/mader/readcounts/ --gene --anno /u/eaudemard/home/data/genome/human/gff3/Homo_sapiens.GRCh38.88.gff3 --subgroup ER --query y
#epcy kal2mat -d ./data/design/mader/all/design.tsv -o ./data/mader/readcounts_bs/ --bs 100 --gene --anno /u/eaudemard/home/data/genome/human/gff3/Homo_sapiens.GRCh38.88.gff3 --subgroup ER --query y

#for design in ${designs}
#do
#  epcy pred_rna -d ./data/design/mader/$design/design.tsv -m ./data/mader/readcounts/readcounts.xls --cpm --log -l 0 -e 0 -o ./data/design/mader/$design/readcounts/ -t 40 -b 10 --kal --min_bw 0.2 --query ${design} --randomseed 42 --full
#  epcy qc -p ./data/design/mader/$design/readcounts/predictive_capability.xls -o ./data/design/mader/$design/readcounts/ --l2fc --min_bw 0.2
#done


epcy profile_rna --query ER -d ./data/design/mader/ER/design.tsv -m ./data/mader/readcounts/readcounts.xls --cpm --log --kal -o ./data/design/mader/ER/readcounts/profile --ids ENSG00000091831 ENSG00000197308 ENSG00000082175 ENSG00000120262 ENSG00000186910 ENSG00000163362 ENSG00000173467 ENSG00000214100 ENSG00000134830 ENSG00000282728 ENSG00000050628
#epcy profile_rna --query HER2 -d ./data/design/mader/HER2/design.tsv -m ./data/mader/readcounts/readcounts.xls --cpm --log --kal -o ./data/design/mader/HER2/readcounts/profile --ids ENSG00000141736 ENSG00000131748 ENSG00000161395 ENSG00000141738 ENSG00000141741 ENSG00000173991
epcy profile_rna --query TNBC -d ./data/design/mader/TNBC/design.tsv -m ./data/mader/readcounts/readcounts.xls --cpm --log --kal -o ./data/design/mader/TNBC/readcounts/profile --ids ENSG00000172425 ENSG00000183888 ENSG00000259459 ENSG00000258910 ENSG00000233078 ENSG00000240800 ENSG00000235687 ENSG00000139865 ENSG00000006071
#epcy profile_rna --query relapse_ER -d ./data/design/mader/relapse_ER/design.tsv -m ./data/mader/readcounts/readcounts.xls --cpm --log --kal -o ./data/design/mader/relapse_ER/readcounts/profile --ids ENSG00000171962
#epcy profile_rna --query relapse_HER2 -d ./data/design/mader/relapse_HER2/design.tsv -m ./data/mader/readcounts/readcounts.xls --cpm --log --kal -o ./data/design/mader/relapse_HER2/readcounts/profile --ids ENSG00000117262 ENSG00000171617 ENSG00000136695
epcy profile_rna --query relapse_TNBC -d ./data/design/mader/relapse_TNBC/design.tsv -m ./data/mader/readcounts/readcounts.xls --cpm --log --kal -o ./data/design/mader/relapse_TNBC/readcounts/profile --ids ENSG00000186205 ENSG00000100865 ENSG00000255750 ENSG00000135925 ENSG00000059573 ENSG00000156587 ENSG00000139192

epcy explore -p ./data/design/mader/ER/readcounts/predictive_capability.xls -s ./data/design/mader/ER/readcounts/subgroup_predicted.xls -o ./data/design/mader/ER/readcounts/explore/ --query ER -d ./data/design/mader/ER/design.tsv




