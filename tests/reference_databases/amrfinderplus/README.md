# Tiny test database preparation

This database was built manually following these steps:

## Step 1 - Download the database
```bash
$ docker run -v .:/data -it quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0 bash
$ amrfinder_update -d amrfinderdb
```

## Step 2 - Reduce size
The database consists of several files. Most are small and kept intact, but the HMM profiles library and proteins file were reduced in size.

### HMM reduction
Using HMMER 3.4 (Aug 2023); http://hmmer.org/
```bash
$ hmmfetch TO_DELETE/AMR.LIB 154989_fosA-NCBIFAM > AMR.LIB
$ hmmpress AMR.LIB
```

### Proteins reduction
The following sequence was selected from the `AMRProt` file:
`>740629697|WP_038415208.1|1|1|fosA|fos_A_A2|thiol_transferase|2|FOSFOMYCIN|FOSFOMYCIN|fosfomycin_resistance_glutathione_transferase_FosA`

A BLAST database was then generated using BLAST version 2.17.0+:
```bash
$ makeblastdb -in AMRProt -dbtype prot
```