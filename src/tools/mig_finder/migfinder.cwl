cwlVersion: v1.2
class: CommandLineTool

label: "MIG_finder runner"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/microbiome-informatics/migfinder"

baseCommand: ["run_migfinder.py"]

inputs:
  query_fasta:
    type: File
    format: edam:format_1929 # FASTA
    inputBinding:
      prefix: "-f"

outputs:
  mig_prediction:
    type: Directory?
    outputBinding:
      glob: data


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
