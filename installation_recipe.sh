git clone https://github.com/mmiladi/MutARNA.git
cd MutARNA
conda env create -f conda-env-MutaRNA.yml -n MutaRNA-env
source activate MutaRNA-env

# Clone R-chie repository
git clone https://github.com/mmiladi/r-chie.git
cd r-chie

# Install R-chie's R4RNA library within the conda env R's library!
my_R_LIB_PATH=$(dirname $(which R))/../lib/R/library/
R -e 'install.packages("R4RNA", repos = NULL, destdir="'$my_R_LIB_PATH'")'

# test run
cd ..
python bin/MutaRNA-plot.py --fasta-wildtype data/sample0.fa --SNP-tag G3C
