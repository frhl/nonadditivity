# Restrict sample IDs to the Europeans
# Run the simulate phenotypes code
# Upload the phenotypes file to the RAP


# to the collection of qced Europeans - the naming of the file 
# is misleading
readonly ancestry="eur" # afr
readonly RAP_location="/wes_ko_ukbb/data/saige/step0/vr_20k/"
readonly sampleIDs="UKB.array.${ancestry}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
readonly sparseGRM="UKB.array.${ancestry}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"

#conda activate random_phenos_v1
readonly prevalence="0.0005,0.001,0.005,0.01,0.1"
readonly heritabilities="0.10"


Rscript 00_make_sim_phenos.R --samples ${sampleIDs} --prevalence ${prevalence} --heritabilities ${heritabilities} --seed 1236 --out 7feb25_uncorrelated_phenos.txt.gz --n_reps=10

