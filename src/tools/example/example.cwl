cwlVersion: v1.2
class: CommandLineTool

label: Example

# Microbiome Informatics
# hints:
#   DockerRequirement:
#     dockerPull: "microbiomeinformatics/mobilomify-<tool>:v<vesrion>"

baseCommand: echo

inputs:
  message:
    type: string
    inputBinding:
      position: 1

arguments: []

stdout: output.txt

outputs:
  example_out:
    type: File
    outputBinding:
      glob: output.txt

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