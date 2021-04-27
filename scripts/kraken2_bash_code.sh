# kraken2 database construction. Seven database were constructed in kraken2.   
# First is to download taxonomic database for all the custom databases.    
  kraken2-build --download-taxonomy --threads 1 --db $DBNAME
# I used the ncbi ITS database as an example for the construction of amplicon database in kraken2.    
  kraken2-build --add-to-library fungi.ITS.fna --db /scratch/cq95/yh7166/db/k2_ncbi_ITS    
  kraken2-build --build --db /scratch/cq95/yh7166/db/k2_ncbi_ITS --threads 14   
  
# Below is a example output for successful construction of the amplicon database.    
  Found 11449/11456 targets, searched through 750580161 accession IDs, search complete.    
  lookup_accession_numbers: 7/11456 accession numbers remain unmapped, see unmapped.txt in DB directory    
  Taxonomy parsed and converted.    
  CHT created with 14 bits reserved for taxid.    
  Completed processing of 11455 sequences, 6742938 bp    
  Writing data to disk... complete.
  
# Here is the general commands for kraken2 analysis of quality filtered nanopore amplicon data   
# I used 16 cpus and 150 Gb memory for the analysis   
  module load kraken2/2.0.8   
  export KRAKEN2_DB_PATH="/scratch/cq95/yh7166/db:"   
  kraken2 --db $db --threads $threads ./barcode01.guppy360.ITS.fasta > ./barcode01_amplicon.${dbname}_output