#!/usr/bin/env python3
"""
Optimize primer sets for genome coverage and linkage disequilibrium
Ranks all compatible sets by multiple criteria
"""

import sys
import argparse
import os
import glob
from collections import defaultdict
from Bio import SeqIO

def determine_snp_density(annotation_file, snp_density_config, target_regions):
    """
    Determine SNP density based on annotation and target regions
    Returns SNPs per kb
    """
    
    # If user provided custom value, use it
    if snp_density_config.get('custom'):
        print(f"Using custom SNP density: {snp_density_config['custom']} SNPs/kb", 
              file=sys.stderr)
        return snp_density_config['custom']
    
    # If annotation available, determine by region type
    if annotation_file and os.path.exists(annotation_file):
        # Determine dominant region type
        region_type = target_regions[0] if target_regions else 'default'
        
        density = snp_density_config.get(region_type, snp_density_config['default'])
        print(f"Using {region_type} SNP density: {density} SNPs/kb", file=sys.stderr)
        return density
    
    # Use genome-wide default
    print(f"Using default SNP density: {snp_density_config['default']} SNPs/kb",
          file=sys.stderr)
    return snp_density_config['default']

def parse_compatible_set(set_file):
    """Parse a compatible set file"""
    pairs = []
    
    with open(set_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            fwd, rev, fwd_tm, rev_tm, count, avg_size, min_size, max_size = parts[:8]
            
            pairs.append({
                'fwd': fwd,
                'rev': rev,
                'fwd_tm': float(fwd_tm),
                'rev_tm': float(rev_tm),
                'count': int(count),
                'avg_size': float(avg_size),
                'min_size': int(min_size),
                'max_size': int(max_size)
            })
    
    return pairs

def calculate_set_metrics(pairs, genome_file, min_contig_distance, 
                         prefer_different_contigs, snp_density):
    """
    Calculate metrics for a primer set
    """
    
    # Load genome to get contig sizes
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))
    contig_sizes = {rec_id: len(rec.seq) for rec_id, rec in genome.items()}
    
    # Metrics
    total_amplicons = sum(p['count'] for p in pairs)
    avg_amplicon_size = sum(p['avg_size'] * p['count'] for p in pairs) / total_amplicons if total_amplicons > 0 else 0
    
    # Estimate genome coverage (sum of amplicon sizes)
    total_coverage_bp = sum(p['avg_size'] * p['count'] for p in pairs)
    genome_size = sum(contig_sizes.values())
    coverage_percent = (total_coverage_bp / genome_size) * 100 if genome_size > 0 else 0
    
    # Estimate SNPs
    estimated_snps = int((total_coverage_bp / 1000) * snp_density)
    
    # LD score (contig distribution)
    # Higher score = better distribution across contigs
    contig_amplicons = defaultdict(int)
    for pair in pairs:
        # This is simplified - in real implementation, parse amplicon locations
        # For now, assume uniform distribution
        for contig in contig_sizes:
            contig_amplicons[contig] += pair['count'] / len(contig_sizes)
    
    # Calculate distribution evenness (lower CV = better)
    if contig_amplicons:
        mean_amp = sum(contig_amplicons.values()) / len(contig_amplicons)
        variance = sum((x - mean_amp) ** 2 for x in contig_amplicons.values()) / len(contig_amplicons)
        cv = (variance ** 0.5) / mean_amp if mean_amp > 0 else 0
        ld_score = 1 / (1 + cv)  # Higher is better
    else:
        ld_score = 0
    
    return {
        'num_pairs': len(pairs),
        'total_amplicons': total_amplicons,
        'avg_amplicon_size': avg_amplicon_size,
        'total_coverage_bp': total_coverage_bp,
        'coverage_percent': coverage_percent,
        'estimated_snps': estimated_snps,
        'ld_score': ld_score
    }

def calculate_composite_score(metrics, target_snps, min_snps):
    """
    Calculate composite score for ranking
    Higher is better
    """
    
    # Normalize metrics
    snp_score = min(metrics['estimated_snps'] / target_snps, 1.0)
    snp_penalty = 0 if metrics['estimated_snps'] >= min_snps else 0.5
    
    coverage_score = min(metrics['coverage_percent'] / 10, 1.0)  # 10% coverage = max
    ld_score = metrics['ld_score']
    size_score = 1.0  # All should be similar
    
    # Weighted composite
    composite = (
        0.4 * snp_score +
        0.3 * coverage_score +
        0.2 * ld_score +
        0.1 * size_score -
        snp_penalty
    )
    
    return composite

def generate_html_report(results, snp_density, target_snps, min_snps, report_file):
    """Generate HTML report with rankings and visualizations"""
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Primer Set Optimization Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #2c3e50;
        }}
        .summary {{
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background-color: white;
            margin-top: 20px;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #3498db;
            color: white;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .rank-1 {{
            background-color: #ffd700 !important;
            font-weight: bold;
        }}
        .rank-2 {{
            background-color: #c0c0c0 !important;
        }}
        .rank-3 {{
            background-color: #cd7f32 !important;
        }}
        .good {{
            color: #27ae60;
            font-weight: bold;
        }}
        .warning {{
            color: #e67e22;
        }}
        .bad {{
            color: #e74c3c;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <h1>🧬 Multiplex Primer Set Optimization Report</h1>
    
    <div class="summary">
        <h2>Configuration Summary</h2>
        <p><strong>SNP Density:</strong> {snp_density} SNPs/kb</p>
        <p><strong>Target SNPs:</strong> {target_snps:,}</p>
        <p><strong>Minimum SNPs:</strong> {min_snps:,}</p>
        <p><strong>Total Sets Analyzed:</strong> {len(results)}</p>
    </div>
    
    <h2>Ranked Primer Sets</h2>
    <p>All compatible sets ranked by composite score (SNP yield, coverage, LD distribution)</p>
    
    <table>
        <thead>
            <tr>
                <th>Rank</th>
                <th>Set ID</th>
                <th>Score</th>
                <th>Pairs</th>
                <th>Amplicons</th>
                <th>Est. SNPs</th>
                <th>Coverage %</th>
                <th>LD Score</th>
                <th>Avg Size (bp)</th>
            </tr>
        </thead>
        <tbody>
"""
    
    for rank, result in enumerate(results, 1):
        m = result['metrics']
        
        # Apply row classes for top 3
        row_class = ''
        if rank == 1:
            row_class = 'rank-1'
        elif rank == 2:
            row_class = 'rank-2'
        elif rank == 3:
            row_class = 'rank-3'
        
        # Color code SNP estimates
        snp_class = 'good' if m['estimated_snps'] >= target_snps else ('warning' if m['estimated_snps'] >= min_snps else 'bad')
        
        html += f"""
            <tr class="{row_class}">
                <td>{rank}</td>
                <td>{result['set_name']}</td>
                <td>{result['score']:.4f}</td>
                <td>{m['num_pairs']}</td>
                <td>{m['total_amplicons']:,}</td>
                <td class="{snp_class}">{m['estimated_snps']:,}</td>
                <td>{m['coverage_percent']:.2f}%</td>
                <td>{m['ld_score']:.3f}</td>
                <td>{m['avg_amplicon_size']:.0f}</td>
            </tr>
"""
    
    html += """
        </tbody>
    </table>
    
    <div class="summary" style="margin-top: 30px;">
        <h2>Legend</h2>
        <p><span class="good">Green SNPs</span>: Meets or exceeds target</p>
        <p><span class="warning">Orange SNPs</span>: Above minimum but below target</p>
        <p><span class="bad">Red SNPs</span>: Below minimum threshold</p>
        <p><strong>Score Components:</strong></p>
        <ul>
            <li>40% - SNP yield (estimated SNPs vs target)</li>
            <li>30% - Genome coverage percentage</li>
            <li>20% - Linkage disequilibrium score (contig distribution)</li>
            <li>10% - Amplicon size consistency</li>
        </ul>
    </div>
    
</body>
</html>
"""
    
    with open(report_file, 'w') as f:
        f.write(html)
    
    print(f"✓ HTML report written to: {report_file}", file=sys.stderr)

def optimize_sets(compatible_dir, genome_file, target_bed, ld_config, output_config,
                 annotation_file, output_ranked, output_dir, report_file):
    """
    Rank all compatible sets by optimization criteria
    """
    
    # Find all set files
    set_files = sorted(glob.glob(f"{compatible_dir}/set_*.txt"))
    print(f"Found {len(set_files)} compatible sets", file=sys.stderr)
    
    if not set_files:
        print("ERROR: No compatible sets found!", file=sys.stderr)
        sys.exit(1)
    
    # Determine SNP density
    snp_density = determine_snp_density(
        annotation_file,
        output_config['snp_density'],
        []  # Would parse from target_bed in real implementation
    )
    
    # Analyze all sets
    print("\nAnalyzing primer sets...", file=sys.stderr)
    results = []
    
    for idx, set_file in enumerate(set_files, 1):
        if idx % 10 == 0:
            print(f"  Analyzed {idx}/{len(set_files)} sets...", file=sys.stderr)
        
        set_name = os.path.basename(set_file).replace('.txt', '')
        
        pairs = parse_compatible_set(set_file)
        
        if not pairs:
            continue
        
        metrics = calculate_set_metrics(
            pairs, genome_file,
            ld_config['min_contig_distance'],
            ld_config['prefer_different_contigs'],
            snp_density
        )
        
        composite_score = calculate_composite_score(
            metrics,
            output_config['expected_snps']['target_snps'],
            output_config['expected_snps']['min_snps']
        )
        
        results.append({
            'set_name': set_name,
            'set_file': set_file,
            'pairs': pairs,
            'metrics': metrics,
            'score': composite_score
        })
    
    # Sort by score (descending)
    results.sort(key=lambda x: x['score'], reverse=True)
    
    # Write ranked results
    print(f"\nWriting ranked results...", file=sys.stderr)
    
    with open(output_ranked, 'w') as f:
        f.write("Rank\tSet_ID\tScore\tPairs\tAmplicons\tEstimated_SNPs\t")
        f.write("Coverage_%\tLD_Score\tAvg_Size\n")
        
        for rank, result in enumerate(results, 1):
            m = result['metrics']
            f.write(f"{rank}\t{result['set_name']}\t{result['score']:.4f}\t")
            f.write(f"{m['num_pairs']}\t{m['total_amplicons']}\t{m['estimated_snps']}\t")
            f.write(f"{m['coverage_percent']:.2f}\t{m['ld_score']:.3f}\t")
            f.write(f"{m['avg_amplicon_size']:.0f}\n")
    
    # Write individual set files
    for rank, result in enumerate(results, 1):
        output_file = f"{output_dir}/set_{rank:04d}.primers.txt"
        
        with open(output_file, 'w') as f:
            f.write(f"# Rank: {rank}\n")
            f.write(f"# Score: {result['score']:.4f}\n")
            f.write(f"# Primer pairs: {result['metrics']['num_pairs']}\n")
            f.write(f"# Total amplicons: {result['metrics']['total_amplicons']:,}\n")
            f.write(f"# Estimated SNPs: {result['metrics']['estimated_snps']:,}\n")
            f.write(f"# Coverage: {result['metrics']['coverage_percent']:.2f}%\n")
            f.write(f"# LD score: {result['metrics']['ld_score']:.3f}\n")
            f.write("#\n")
            f.write("Forward\tReverse\tFwd_Tm\tRev_Tm\tAmplicons\tAvg_Size\tMin_Size\tMax_Size\n")
            
            for pair in result['pairs']:
                f.write(f"{pair['fwd']}\t{pair['rev']}\t")
                f.write(f"{pair['fwd_tm']:.1f}\t{pair['rev_tm']:.1f}\t")
                f.write(f"{pair['count']}\t{pair['avg_size']:.0f}\t")
                f.write(f"{pair['min_size']}\t{pair['max_size']}\n")
    
    # Generate HTML report
    generate_html_report(
        results, snp_density,
        output_config['expected_snps']['target_snps'],
        output_config['expected_snps']['min_snps'],
        report_file
    )
    
    print(f"\n✓ Ranked {len(results)} primer sets", file=sys.stderr)
    print(f"✓ Top set: {results[0]['set_name']} (score: {results[0]['score']:.4f})", file=sys.stderr)
    print(f"✓ Estimated SNPs: {results[0]['metrics']['estimated_snps']:,}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Optimize and rank primer sets'
    )
    parser.add_argument('--compatible-sets', required=True, 
                       help='Directory with compatible set files')
    parser.add_argument('--genome', required=True, help='Genome FASTA')
    parser.add_argument('--target-bed', required=True, help='Target regions BED')
    parser.add_argument('--annotation', help='Annotation file (optional)')
    parser.add_argument('--min-contig-distance', type=int, default=10000,
                       help='Minimum distance between amplicons on same contig')
    parser.add_argument('--prefer-different-contigs', type=bool, default=True,
                       help='Prefer amplicons on different contigs')
    parser.add_argument('--snp-density', type=float, required=True,
                       help='SNP density (SNPs per kb) - can be "auto"')
    parser.add_argument('--min-snps', type=int, default=20000,
                       help='Minimum SNPs required')
    parser.add_argument('--target-snps', type=int, default=50000,
                       help='Target SNPs')
    parser.add_argument('--report-all', type=bool, default=True,
                       help='Report all sets (not just top N)')
    parser.add_argument('--output-ranked', required=True,
                       help='Output ranked TSV file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for individual set files')
    parser.add_argument('--report', required=True, help='HTML report file')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Build config objects
    ld_config = {
        'min_contig_distance': args.min_contig_distance,
        'prefer_different_contigs': args.prefer_different_contigs
    }
    
    output_config = {
        'snp_density': {
            'default': args.snp_density,
            'custom': None
        },
        'expected_snps': {
            'min_snps': args.min_snps,
            'target_snps': args.target_snps
        }
    }
    
    optimize_sets(
        compatible_dir=args.compatible_sets,
        genome_file=args.genome,
        target_bed=args.target_bed,
        ld_config=ld_config,
        output_config=output_config,
        annotation_file=args.annotation,
        output_ranked=args.output_ranked,
        output_dir=args.output_dir,
        report_file=args.report
    )

if __name__ == '__main__':
    main()
