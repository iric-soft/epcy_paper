
#To build data/other
#Rscript --vanilla ./src/script/other/biotype.r
#Rscript --vanilla ./src/script/other/create_gene_length.r

# read count to tpm
#Rscript --vanilla ./src/script/other/create_matrix_tpm.r

#To build data/TCGA_BRCA
#Rscript --vanilla ./src/script/other/brca_matrix_readcounts.r


design="62_CK"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="9_EVI1"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="62_Inter"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="inv16"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="33_MLL"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="13_Mono5"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="126_Normal"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="5_NUP98NSD1"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="30_t15_17"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="18_t8_21"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Query/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Ref/g' ./data/design/leucegene/${design}/design.tsv

design="13_Tri8"
cp ./data/design/leucegene/all/design.tsv ./data/design/leucegene/${design}/
sed -i -e 's/CK/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/EVI1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Inter/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/inv16/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/MLL/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Mono5/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Normal/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/NUP98NSD1/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t15_17/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/t8_21/Ref/g' ./data/design/leucegene/${design}/design.tsv
sed -i -e 's/Tri8/Query/g' ./data/design/leucegene/${design}/design.tsv
