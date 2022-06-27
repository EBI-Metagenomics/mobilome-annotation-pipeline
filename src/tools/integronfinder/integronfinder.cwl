cwlVersion: v1.2
class: CommandLineTool

label: "IntegronFinder runner"


requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/microbiome-informatics/integronfinder"

baseCommand: ["integron_finder"]

arguments: [--union-integrases, --mute, --local-max, --cpu, '4', --func-annot, --outdir, 'output']

inputs:
  replicon:
    type: File
    format: edam:format_1929 # FASTA
    inputBinding:
      position: 9

outputs:
  out_path:
    type: Directory?
    outputBinding:
      glob: 'output'


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
