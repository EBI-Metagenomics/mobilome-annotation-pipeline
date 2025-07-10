process GENERATE_GFF3_OUTPUT {
    tag "$meta.id"
    
    publishDir "${params.outdir}/final_gff3", mode: 'copy'
    
    input:
    tuple val(meta), path(quality_assessed_ice)
    tuple val(meta), path(integrated_ice_info)
    path genome_fasta
    
    output:
    tuple val(meta), path("${prefix}_ice_elements.gff3"), emit: ice_gff3
    tuple val(meta), path("${prefix}_ice_features.gff3"), emit: ice_features_gff3
    tuple val(meta), path("${prefix}_ice_summary.tsv"), emit: ice_summary
    path "versions.yml", emit: versions
    
    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/usr/bin/env python3
    
    import json
    import pandas as pd
    from Bio import SeqIO
    import sys
    from datetime import datetime
    
    def generate_gff3_header():
        '''Generate GFF3 header with metadata'''
        
        header_lines = [
            "##gff-version 3",
            f"##date {datetime.now().strftime('%Y-%m-%d')}",
            "##source icefinder2-nextflow",
            "##feature-ontology http://song.cvs.sourceforge.net/viewvc/song/ontology/sofa.obo?revision=1.269",
            "##species Unknown",
            ""
        ]
        
        return header_lines
    
    def generate_ice_elements_gff3(quality_assessed_ice, integrated_info):
        '''Generate GFF3 for ICE elements'''
        
        gff3_lines = generate_gff3_header()
        
        for _, ice in quality_assessed_ice.iterrows():
            if ice['assessment_status'] != 'PASS':
                continue
                
            system_id = ice['system_id']
            contig = ice['contig']
            start = ice['start']
            end = ice['end']
            ice_type = ice['ice_type']
            quality_score = ice['quality_score']
            
            # Get additional information from integrated data
            ice_info = integrated_info.get(system_id, {})
            basic_info = ice_info.get('basic_info', {})
            orit_info = ice_info.get('orit_info', {})
            trna_info = ice_info.get('trna_info', {})
            dr_info = ice_info.get('direct_repeats', {})
            
            # Create ICE element entry
            attributes = [
                f"ID={system_id}",
                f"Name={system_id}",
                f"ice_type={ice_type}",
                f"length={ice['length']}",
                f"quality_score={quality_score}",
                f"confidence={ice['confidence']}",
                f"gc_content={basic_info.get('gc_content', 0)}",
                f"num_genes={ice['num_integrase'] + ice['num_relaxase'] + ice['num_t4ss'] + ice['num_t4cp'] + ice['num_replication'] + ice['num_transfer']}",
                f"has_integrase={ice['num_integrase'] > 0}",
                f"has_relaxase={ice['num_relaxase'] > 0}",
                f"has_t4ss={ice['num_t4ss'] > 0}",
                f"has_orit={ice['has_orit']}",
                f"num_trnas={ice['num_trnas']}",
                f"mobilization_types={ice['mobilization_types']}",
                f"t4ss_types={ice['t4ss_types']}"
            ]
            
            if orit_info:
                attributes.extend([
                    f"orit_identity={orit_info.get('identity_percent', 0)}",
                    f"orit_coverage={orit_info.get('coverage', 0)}",
                    f"orit_start={orit_info.get('query_start', 0)}",
                    f"orit_end={orit_info.get('query_end', 0)}"
                ])
            
            if trna_info.get('closest_trna'):
                closest_trna = trna_info['closest_trna']
                attributes.append(f"closest_trna={closest_trna.get('type', 'unknown')}")
            
            if dr_info.get('boundaries'):
                dr_boundaries = dr_info['boundaries']
                attributes.extend([
                    f"attL={dr_boundaries.get('attL_start', 0)}-{dr_boundaries.get('attL_end', 0)}",
                    f"attR={dr_boundaries.get('attR_start', 0)}-{dr_boundaries.get('attR_end', 0)}",
                    f"dr_length={dr_boundaries.get('repeat_length', 0)}"
                ])
            
            # Create GFF3 line
            gff3_line = f"{contig}\\ticefinder2\\tintegrative_conjugative_element\\t{start}\\t{end}\\t{quality_score}\\t.\\t.\\t{';'.join(attributes)}"
            gff3_lines.append(gff3_line)
        
        return gff3_lines
    
    def generate_ice_features_gff3(quality_assessed_ice, integrated_info):
        '''Generate detailed GFF3 for ICE features'''
        
        gff3_lines = generate_gff3_header()
        
        for _, ice in quality_assessed_ice.iterrows():
            if ice['assessment_status'] != 'PASS':
                continue
                
            system_id = ice['system_id']
            contig = ice['contig']
            
            # Get feature information
            ice_info = integrated_info.get(system_id, {})
            features = ice_info.get('features', [])
            orit_info = ice_info.get('orit_info', {})
            trna_info = ice_info.get('trna_info', {})
            dr_info = ice_info.get('direct_repeats', {})
            
            # Add ICE element as parent
            ice_attributes = f"ID={system_id};Name={system_id};ice_type={ice['ice_type']}"
            ice_line = f"{contig}\\ticefinder2\\tintegrative_conjugative_element\\t{ice['start']}\\t{ice['end']}\\t{ice['quality_score']}\\t.\\t.\\t{ice_attributes}"
            gff3_lines.append(ice_line)
            
            # Add individual features
            for i, feature in enumerate(features):
                feature_id = f"{system_id}_feature_{i+1}"
                gene_name = feature['gene_name']
                start = feature['start']
                end = feature['end']
                strand = feature['strand']
                
                # Determine feature type based on gene name
                if any(pattern in gene_name for pattern in ['integrase', 'Phage_integrase', 'rve']):
                    feature_type = 'integrase'
                elif any(pattern in gene_name for pattern in ['Relaxase', 'MOB']):
                    feature_type = 'relaxase'
                elif any(pattern in gene_name for pattern in ['T4SS', 'FATA', 'FA_']):
                    feature_type = 'type_IV_secretion_system'
                elif 't4cp' in gene_name.lower():
                    feature_type = 'type_IV_coupling_protein'
                elif any(pattern in gene_name for pattern in ['Rep', 'Prim-Pol']):
                    feature_type = 'replication_protein'
                else:
                    feature_type = 'CDS'
                
                feature_attributes = [
                    f"ID={feature_id}",
                    f"Parent={system_id}",
                    f"Name={gene_name}",
                    f"gene_name={gene_name}",
                    f"feature_type={feature_type}"
                ]
                
                feature_line = f"{contig}\\tmacsyfinder\\t{feature_type}\\t{start}\\t{end}\\t.\\t{strand}\\t.\\t{';'.join(feature_attributes)}"
                gff3_lines.append(feature_line)
            
            # Add oriT if present
            if orit_info:
                orit_id = f"{system_id}_oriT"
                orit_start = ice['start'] + orit_info.get('query_start', 1) - 1
                orit_end = ice['start'] + orit_info.get('query_end', 1) - 1
                
                orit_attributes = [
                    f"ID={orit_id}",
                    f"Parent={system_id}",
                    f"Name=oriT",
                    f"identity={orit_info.get('identity_percent', 0)}",
                    f"coverage={orit_info.get('coverage', 0)}",
                    f"evalue={orit_info.get('evalue', 0)}",
                    f"reference={orit_info.get('orit_id', 'unknown')}"
                ]
                
                orit_line = f"{contig}\\ticefinder2\\torigin_of_transfer\\t{orit_start}\\t{orit_end}\\t{orit_info.get('h_value', 0)}\\t.\\t.\\t{';'.join(orit_attributes)}"
                gff3_lines.append(orit_line)
            
            # Add tRNAs
            for j, trna in enumerate(trna_info.get('associated_trnas', [])):
                trna_id = f"{system_id}_tRNA_{j+1}"
                trna_attributes = [
                    f"ID={trna_id}",
                    f"Parent={system_id}",
                    f"Name={trna['type']}",
                    f"anticodon={trna.get('anticodon', 'unknown')}"
                ]
                
                trna_line = f"{contig}\\ttRNAscan-SE\\ttRNA\\t{trna['start']}\\t{trna['end']}\\t.\\t{trna['strand']}\\t.\\t{';'.join(trna_attributes)}"
                gff3_lines.append(trna_line)
            
            # Add direct repeats (attL and attR)
            if dr_info.get('boundaries'):
                dr_boundaries = dr_info['boundaries']
                
                # attL (left attachment site)
                attL_id = f"{system_id}_attL"
                attL_attributes = [
                    f"ID={attL_id}",
                    f"Parent={system_id}",
                    f"Name=attL",
                    f"repeat_length={dr_boundaries.get('repeat_length', 0)}"
                ]
                attL_line = f"{contig}\\ticefinder2\\tattachment_site\\t{dr_boundaries.get('attL_start', 0)}\\t{dr_boundaries.get('attL_end', 0)}\\t.\\t.\\t.\\t{';'.join(attL_attributes)}"
                gff3_lines.append(attL_line)
                
                # attR (right attachment site)
                attR_id = f"{system_id}_attR"
                attR_attributes = [
                    f"ID={attR_id}",
                    f"Parent={system_id}",
                    f"Name=attR",
                    f"repeat_length={dr_boundaries.get('repeat_length', 0)}"
                ]
                attR_line = f"{contig}\\ticefinder2\\tattachment_site\\t{dr_boundaries.get('attR_start', 0)}\\t{dr_boundaries.get('attR_end', 0)}\\t.\\t.\\t.\\t{';'.join(attR_attributes)}"
                gff3_lines.append(attR_line)
        
        return gff3_lines
    
    def generate_ice_summary(quality_assessed_ice):
        '''Generate ICE summary statistics'''
        
        summary_data = []
        
        # Overall statistics
        total_ice = len(quality_assessed_ice)
        passed_ice = len(quality_assessed_ice[quality_assessed_ice['assessment_status'] == 'PASS'])
        
        # Type distribution
        type_counts = quality_assessed_ice['ice_type'].value_counts().to_dict()
        
        # Size statistics
        sizes = quality_assessed_ice['length']
        size_stats = {
            'min_size': sizes.min(),
            'max_size': sizes.max(),
            'mean_size': round(sizes.mean()),
            'median_size': round(sizes.median())
        }
        
        # Feature statistics
        feature_stats = {
            'with_integrase': len(quality_assessed_ice[quality_assessed_ice['num_integrase'] > 0]),
            'with_relaxase': len(quality_assessed_ice[quality_assessed_ice['num_relaxase'] > 0]),
            'with_t4ss': len(quality_assessed_ice[quality_assessed_ice['num_t4ss'] > 0]),
            'with_orit': len(quality_assessed_ice[quality_assessed_ice['has_orit'] == True]),
            'with_trnas': len(quality_assessed_ice[quality_assessed_ice['num_trnas'] > 0])
        }
        
        # Quality distribution
        quality_bins = pd.cut(quality_assessed_ice['quality_score'], 
                            bins=[0, 0.5, 0.8, 1.0], 
                            labels=['low', 'medium', 'high'])
        quality_distribution = quality_bins.value_counts().to_dict()
        
        summary_data.append({
            'metric': 'total_ice_elements',
            'value': total_ice,
            'description': 'Total number of ICE elements processed'
        })
        
        summary_data.append({
            'metric': 'passed_quality_assessment',
            'value': passed_ice,
            'description': 'Number of ICE elements that passed quality assessment'
        })
        
        summary_data.append({
            'metric': 'pass_rate',
            'value': f"{passed_ice/total_ice*100:.1f}%" if total_ice > 0 else "0%",
            'description': 'Percentage of ICE elements that passed quality assessment'
        })
        
        # Add type distribution
        for ice_type, count in type_counts.items():
            summary_data.append({
                'metric': f'type_{ice_type.lower().replace(" ", "_").replace("-", "_")}',
                'value': count,
                'description': f'Number of {ice_type} elements'
            })
        

        # Add size statistics
        for stat_name, stat_value in size_stats.items():
            summary_data.append({
                'metric': f'size_{stat_name}',
                'value': stat_value,
                'description': f'ICE size {stat_name.replace("_", " ")}'
            })
        
        # Add feature statistics
        for feature_name, feature_count in feature_stats.items():
            summary_data.append({
                'metric': f'elements_{feature_name}',
                'value': feature_count,
                'description': f'Number of ICE elements {feature_name.replace("_", " ")}'
            })
        
        # Add quality distribution
        for quality_level, count in quality_distribution.items():
            summary_data.append({
                'metric': f'quality_{quality_level}',
                'value': count,
                'description': f'Number of ICE elements with {quality_level} quality score'
            })
        
        return summary_data
    
    def main():
        '''Main function for GFF3 output generation'''
        
        # Read input files
        quality_assessed_ice = pd.read_csv('${quality_assessed_ice}', sep='\\t')
        
        with open('${integrated_ice_info}', 'r') as f:
            integrated_info = json.load(f)
        
        print(f"Generating GFF3 output for {len(quality_assessed_ice)} ICE elements")
        
        # Generate ICE elements GFF3
        ice_elements_gff3 = generate_ice_elements_gff3(quality_assessed_ice, integrated_info)
        
        with open('${prefix}_ice_elements.gff3', 'w') as f:
            f.write('\\n'.join(ice_elements_gff3))
        
        # Generate detailed ICE features GFF3
        ice_features_gff3 = generate_ice_features_gff3(quality_assessed_ice, integrated_info)
        
        with open('${prefix}_ice_features.gff3', 'w') as f:
            f.write('\\n'.join(ice_features_gff3))
        
        # Generate ICE summary
        summary_data = generate_ice_summary(quality_assessed_ice)
        
        if summary_data:
            df_summary = pd.DataFrame(summary_data)
            df_summary.to_csv('${prefix}_ice_summary.tsv', sep='\\t', index=False)
        else:
            # Create empty summary file
            with open('${prefix}_ice_summary.tsv', 'w') as f:
                f.write('metric\\tvalue\\tdescription\\n')
        
        print(f"Generated GFF3 files:")
        print(f"  ICE elements: ${prefix}_ice_elements.gff3")
        print(f"  ICE features: ${prefix}_ice_features.gff3")
        print(f"  Summary: ${prefix}_ice_summary.tsv")
        
        # Print summary statistics
        passed_elements = quality_assessed_ice[quality_assessed_ice['assessment_status'] == 'PASS']
        print(f"\\nSummary statistics:")
        print(f"  Total ICE elements: {len(quality_assessed_ice)}")
        print(f"  Passed quality assessment: {len(passed_elements)}")
        print(f"  ICE types: {quality_assessed_ice['ice_type'].value_counts().to_dict()}")
    
    # Run GFF3 generation
    main()
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python --version | sed 's/Python //')"
        pandas: "\$(python -c 'import pandas; print(pandas.__version__)')"
        biopython: "\$(python -c 'import Bio; print(Bio.__version__)')"
    END_VERSIONS
    """
}

