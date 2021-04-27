# application of minimap2 to gold standard alignment database (see fungal_species_consensus_seqs_44.txt) and NCBI database

def minimapmapping(fasta_fn, ref_fn, out_fn, threads):
    command = F"minimap2 -x map-ont -t {threads} {ref_fn} {fasta_fn} -o {out_fn}"
    out = subprocess.getstatusoutput(command)
    
    
# fasta_fn == fasta file for input reads
# ref_fn == database file
# out_fn == output file