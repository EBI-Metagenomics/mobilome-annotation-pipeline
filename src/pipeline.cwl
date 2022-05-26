cwlVersion: v1.2
class: Workflow

label: Mobilomify 

inputs:
  message:
    type: string

steps:
  step1_example:
    label: example
    run: ./tools/example/example.cwl
    in:
      message: message
    out:
      - example_out

outputs:
    example_output:
      outputSource: step1_example/example_out
      type: File

doc: |
    Description of the pipeline

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