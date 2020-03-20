
# 3 High mcc without bs, 3 High mcc with bs, 3 High mcc with all cohort, 3 bagging help
ids="ENSG00000183087 ENSG00000262456 ENSG00000272573 ENSG00000020633 ENSG00000138821 ENSG00000169398 ENSG00000133392 ENSG00000147488 ENSG00000127328 ENSG00000153814 ENSG00000132475 ENSG00000211943"
gff="./data/other/Homo_sapiens.GRCh38.84.gff3"
path_5v5="./data/design/leucegene/5_inv16_vs_5/"
path_all="./data/design/leucegene/28_inv16/"

epcy profile_rna --kal --cpm --log --gene --anno $gff \
                 -d ${path_5v5}/design.tsv \
                 -o ${path_5v5}/profile/ \
                 --ids $ids

epcy profile_rna --kal --cpm --log --gene --anno $gff \
                 -d ${path_5v5}/design.tsv \
                 -o ${path_5v5}/profile_bs10/ \
                 --bs 10 --ids $ids

epcy profile_rna --kal --cpm --log --gene --anno $gff \
                 -d ${path_5v5}/design.tsv \
                 -o ${path_5v5}/profile_bs100/ \
                 --bs 100 --ids $ids

epcy profile_rna --kal --cpm --log --gene --anno $gff \
                -d ${path_all}/design.tsv \
                -o ${path_5v5}/profile_all/ \
                --ids $ids
