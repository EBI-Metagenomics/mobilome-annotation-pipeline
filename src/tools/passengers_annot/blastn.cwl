cwlVersion: v1.2
class: CommandLineTool

label: "blastn vs resfinder"


requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/microbiome-informatics/blast_simple:latest"


baseCommand: ["blastn"]

arguments: [-outfmt, '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident positive qcovs',
          -num_threads, '4',
          -max_target_seqs, '1', 
          -evalue, '0.01']

inputs:
  query_seq:
    type: File
    format: edam:format_1929 # FASTA
    inputBinding:
      prefix: "-query"

  database:
    type: File
    format: http://edamontology.org/operation_2421 # blast database
    inputBinding:
      prefix: "-db"
    secondaryFiles:
      - .ndb
      - .nhr
      - .nin
      - .not
      - .nsq
      - .ntf
      - .nto

stdout: resfinder_blast.out

outputs:
  alignment:
    type: stdout


$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/
$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-http.rdf

s:license: "https://www.apache.org/licenses/LICENSE-2.0"
s:copyrightHolder:
  - class: s:Organization
    s:name: "EMBL - European Bioinformatics Institute"
    s:url: "https://www.ebi.ac.uk"
