process BOUNDARY_DETECTION {
    conda "bioconda::vmatch=2.3.1 bioconda::trnascan-se=2.0.9 bioconda::biopython=1.79"
    
    input:
    tuple val(meta), path(macsyfinder_results)
    tuple val(meta), path(gbk_file)

    output:
    tuple val(meta), path("${meta.id}_ice_boundaries.json"), emit: ice_with_boundaries
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env python3
    
    import json
    import subprocess
    import sys
    from Bio import SeqIO
    import tempfile
    import os
    
    def find_direct_repeats(sequence, min_length=10, max_distance=50000):
        \"\"\"Find direct repeats using vmatch\"\"\"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">sequence\\n{sequence}")
            fasta_file = f.name
        
        try:
            # Build vmatch index
            subprocess.run(["mkvtree", "-db", fasta_file, "-dna", "-pl"], 
                         check=True, capture_output=True)
            
            # Find repeats
            result = subprocess.run([
                "vmatch", "-l", str(min_length), "-best", "10",
                f"{fasta_file}.prj"
            ], capture_output=True, text=True)
            
            repeats = []
            for line in result.stdout.split('\\n'):
                if line.strip() and not line.startswith('#'):
                    fields = line.split()
                    if len(fields) >= 4:
                        length = int(fields[0])
                        pos1 = int(fields[1])
                        pos2 = int(fields[2])
                        if abs(pos2 - pos1) <= max_distance:
                            repeats.append({
                                'length': length,
                                'pos1': pos1,
                                'pos2': pos2,
                                'distance': abs(pos2 - pos1)
                            })
            
            return repeats
            
        finally:
            # Cleanup
            for ext in ['.prj', '.tis', '.ois', '.suf', '.lcp', '.llv', '.des']:
                try:
                    os.unlink(f"{fasta_file},{ext}")
                except:
                    pass
            os.unlink(fasta_file)
    
    def find_trna_genes(sequence):
        \"\"\"Find tRNA genes using tRNAscan-SE\"\"\"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">sequence\\n{sequence}")
            fasta_file = f.name
        
        try:
            result = subprocess.run([
                "tRNAscan-SE", "-B", "-o", "/dev/stdout", fasta_file
            ], capture_output=True, text=True)
            
            trnas = []
            for line in result.stdout.split('\\n'):
                if line.strip() and not line.startswith('Name'):
                    fields = line.split()
                    if len(fields) >= 4:
                        start = int(fields[2])
                        end = int(fields[3])
                        trnas.append({
                            'start': min(start, end),
                            'end': max(start, end),
                            'type': fields[4] if len(fields) > 4 else 'unknown'
                        })
            
            return trnas
            
        finally:
            os.unlink(fasta_file)
    
    def refine_ice_boundaries(ice_data, gbk_file):
        \"\"\"Refine ICE boundaries using direct repeats and tRNA analysis\"\"\"
        # Load GenBank file
        records = list(SeqIO.parse(gbk_file, "genbank"))
        
        refined_ices = []
        
        # Process MacSyFinder results
        results_file = os.path.join("${macsyfinder_results}", "results.json")
        if os.path.exists(results_file):
            with open(results_file, 'r') as f:
                macsyfinder_data = json.load(f)
            
            for system in macsyfinder_data.get('systems', []):
                # Find corresponding sequence
                sequence = None
                for record in records:
                    if record.id in system.get('replicon', ''):
                        sequence = str(record.seq)
                        break
                
                if sequence:
                    # Get ICE region coordinates
                    start = system.get('start', 0)
                    end = system.get('end', len(sequence))
                    
                    # Find direct repeats near boundaries
                    region_start = max(0, start - 5000)
                    region_end = min(len(sequence), end + 5000)
                    region_seq = sequence[region_start:region_end]
                    
                    repeats = find_direct_repeats(region_seq)
                    trnas = find_trna_genes(region_seq)
                    
                    # Refine boundaries based on repeats and tRNAs
                    refined_start = start
                    refined_end = end
                    
                    # Look for direct repeats that could be ICE boundaries
                    for repeat in repeats:
                        if repeat['length'] >= 15:  # Minimum DR length
                            dr_start = region_start + repeat['pos1']
                            dr_end = region_start + repeat['pos2']
                            
                            if abs(dr_start - start) < 2000:
                                refined_start = dr_start
                            if abs(dr_end - end) < 2000:
                                refined_end = dr_end
                    
                    # Check for tRNA integration sites
                    integration_site = None
                    for trna in trnas:
                        trna_pos = region_start + trna['start']
                        if abs(trna_pos - refined_start) < 1000 or abs(trna_pos - refined_end) < 1000:
                            integration_site = {
                                'position': trna_pos,
                                'type': trna['type']
                            }
                            break
                    
                    refined_ice = {
                        'id': f"ICE_{len(refined_ices) + 1}",
                        'replicon': system.get('replicon'),
                        'original_start': start,
                        'original_end': end,
                        'refined_start': refined_start,
                        'refined_end': refined_end,
                        'length': refined_end - refined_start,
                        'direct_repeats': repeats,
                        'integration_site': integration_site,
                        'genes': system.get('genes', []),
                        'system_type': system.get('type', 'ICE')
                    }
                    
                    refined_ices.append(refined_ice)
        
        return refined_ices
    
    # Main execution
    refined_ices = refine_ice_boundaries("${macsyfinder_results}", "${gbk_file}")
    
    # Save results
    with open("${meta.id}_ice_boundaries.json", "w") as f:
        json.dump({
            'sample_id': "${meta.id}",
            'ices': refined_ices,
            'total_ices': len(refined_ices)
        }, f, indent=2)
    
    print(f"Refined boundaries for {len(refined_ices)} ICE elements")
    
    # Write versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    vmatch: "2.3.1"\\n')
        f.write('    trnascan-se: "2.0.9"\\n')
        f.write('    biopython: "1.79"\\n')
    """
}
