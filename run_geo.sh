#module add python/3.6.2

#epcy kal2mat -d ./data/design/GEO/all/design.tsv -o ./data/GEO_PRJNA299579/readcounts_bs/ --bs 20 --gene --anno /u/eaudemard/home/data/genome/human/gff3/Homo_sapiens.GRCh38.88.gff3

designs="LSC_vs_blast LSC_vs_GMP LSC_vs_LMPP LSC_vs_pHSC"

#for design in ${designs}
#do
#  epcy pred_rna -d ./data/design/GEO/$design/design.tsv -m ./data/GEO_PRJNA299579/readcounts_bs/readcounts.xls --cpm --log -l 0 -e 0 -o ./data/design/GEO/$design/readcounts_bs_bagging/ -t 50 -b 10 --bs 20 --kal --min_bw 0.2 --query LSC --randomseed 42
#  epcy qc -p ./data/design/GEO/$design/readcounts_bs_bagging/predictive_capability.xls -o ./data/design/GEO/$design/readcounts_bs_bagging/ --l2fc --min_bw 0.2
#done

#epcy pred_rna --gene -d ./data/design/GEO/all/design.tsv --cpm --log -l 0 -e 0 -o ./data/design/GEO/all/LSC/readcounts/ -t 30 --kal --min_bw 0.2 --query LSC --anno /u/eaudemard/home/data/genome/human/gff3/Homo_sapiens.GRCh38.88.gff3 --randomseed 42
epcy qc -p ./data/design/GEO/all/LSC/readcounts/predictive_capability.xls -o ./data/design/GEO/all/LSC/readcounts/ --l2fc --min_bw 0.2



#epcy profile_rna --query LSC --min_bw 0.2 --bs 20 -d ./data/design/GEO/LSC_vs_blast/design.tsv -m ./data/GEO_PRJNA299579/readcounts_bs/readcounts.xls --cpm --log --kal -o ./data/design/GEO/LSC_vs_blast/readcounts_bs_bagging/profil_LSC --ids ENSG00000121851 ENSG00000133112 ENSG00000114209 ENSG00000161970
#epcy profile_rna --query LSC --min_bw 0.2 --bs 20 -d ./data/design/GEO/LSC_vs_blast/design.tsv -m ./data/GEO_PRJNA299579/readcounts_bs/readcounts.xls --cpm --log --kal -o ./data/design/GEO/LSC_vs_blast/readcounts_bs_bagging/profil_blast --ids ENSG00000181222 ENSG00000142657 ENSG00000149577 ENSG00000105723
#epcy profile_rna --query LSC --min_bw 0.2 --bs 20 -d ./data/design/GEO/LSC_vs_blast/design.tsv -m ./data/GEO_PRJNA299579/readcounts_bs/readcounts.xls --cpm --log --kal -o ./data/design/GEO/LSC_vs_blast/readcounts_bs_bagging/profil_other --ids ENSG00000120217

#epcy profile_rna --query LSC --min_bw 0.2 --bs 20 -d ./data/design/GEO/LSC_vs_pHSC/design.tsv -m ./data/GEO_PRJNA299579/readcounts_bs/readcounts.xls --cpm --log --kal -o ./data/design/GEO/LSC_vs_pHSC/readcounts_bs_bagging/profil_pHSC --ids ENSG00000117400 ENSG00000162367 ENSG00000162599 ENSG00000179639 ENSG00000185630 ENSG00000165821 ENSG00000185650 ENSG00000080603 ENSG00000133026 ENSG00000121104 ENSG00000115641 ENSG00000100234 ENSG00000154783 ENSG00000151090 ENSG00000120279 ENSG00000008311 ENSG00000164946 ENSG00000165152 ENSG00000130723
#epcy profile_rna --query LSC --min_bw 0.2 --bs 20 -d ./data/design/GEO/LSC_vs_pHSC/design.tsv -m ./data/GEO_PRJNA299579/readcounts_bs/readcounts.xls --cpm --log --kal -o ./data/design/GEO/LSC_vs_pHSC/readcounts_bs_bagging/profil_LSC --ids ENSG00000148700

