cwlVersion: v1.2
class: CommandLineTool

label: "blastn vs mge_records database"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/microbiome-informatics/blast_simple:latest"

baseCommand: ["blastn"]

arguments:
  - prefix: -outfmt
    position: 1
    valueFrom: '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident positive qcovs'
  - prefix: -num_threads
    valueFrom: '4'
    position: 2
  - prefix: -max_target_seqs
    valueFrom: '1'
    position: 3
  - prefix: -evalue
    valueFrom: '1E-5'
    position: 4
  - prefix: -dust
    valueFrom: 'no'
    position: 5

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

stdout: mgerecords_blast.out

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
