#!/usr/bin/env python3
"""
Validate primer sets through in silico PCR simulation
Generate visualization reports
"""

import sys
import argparse
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def parse_ranked_sets(ranked_file, top_n=10):
    """Parse ranked sets file and return top N set IDs"""
    sets = []
    
    with open(ranked_file) as f:
        header = f.readline()
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            rank = int(parts[0])
            set_id = parts[1]
            
            if rank <= top_n:
                sets.append(set_id)
    
    return sets

def simulate_pcr(primer_fwd, primer_rev, genome, max_size=500):
    """
    Simulate PCR amplification for a primer pair
    Returns list of amplicons
    """
    amplicons = []
    
    for contig_id, record in genome.items():
        seq = str(record.seq).upper()
        
        # Find forward primer binding sites
        fwd_sites = []
        pos = 0
        while True:
            pos = seq.find(primer_fwd, pos)
            if pos == -1:
                break
            fwd_sites.append(pos)
            pos += 1
        
        # Find reverse primer binding sites (reverse complement)
        rev_rc = str(Seq(primer_rev).reverse_complement())
        rev_sites = []
        pos = 0
        while True:
            pos = seq.find(rev_rc, pos)
            if pos == -1:
                break
            rev_sites.append(pos)
            pos += 1
        
        # Find amplicons (fwd before rev, within size limit)
        for fwd_pos in fwd_sites:
            for rev_pos in rev_sites:
                if fwd_pos < rev_pos:
                    size = rev_pos - fwd_pos + len(primer_rev)
                    
                    if size <= max_size:
                        amplicons.append({
                            'contig': contig_id,
                            'start': fwd_pos,
                            'end': rev_pos + len(primer_rev),
                            'size': size,
                            'sequence': seq[fwd_pos:rev_pos + len(primer_rev)]
                        })
    
    return amplicons

def validate_primer_set(set_file, genome, extension_time, polymerase, output_dir):
    """
    Validate a single primer set through in silico PCR
    """
    
    # Calculate max amplicon size from extension time
    if polymerase.upper() == 'Q5':
        max_size = int(extension_time / 30 * 1000)
    else:
        max_size = int(extension_time / 60 * 1000)
    
    # Parse primer set
    primers = []
    with open(set_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                primers.append({
                    'fwd': parts[0],
                    'rev': parts[1]
                })
    
    # Simulate PCR for each pair
    all_amplicons = []
    
    for i, pair in enumerate(primers):
        amplicons = simulate_pcr(pair['fwd'], pair['rev'], genome, max_size)
        
        for amp in amplicons:
            amp['pair_id'] = i
            amp['fwd'] = pair['fwd']
            amp['rev'] = pair['rev']
        
        all_amplicons.extend(amplicons)
    
    return all_amplicons, primers

def generate_bed_file(amplicons, output_file):
    """Generate BED file of simulated amplicons"""
    
    with open(output_file, 'w') as f:
        for i, amp in enumerate(amplicons):
            f.write(f"{amp['contig']}\t{amp['start']}\t{amp['end']}\t")
            f.write(f"amplicon_{i+1}_pair_{amp['pair_id']+1}\t0\t+\n")

def generate_fasta_file(amplicons, output_file):
    """Generate FASTA file of simulated amplicon sequences"""
    
    with open(output_file, 'w') as f:
        for i, amp in enumerate(amplicons):
            f.write(f">amplicon_{i+1}_pair_{amp['pair_id']+1}_{amp['contig']}:{amp['start']}-{amp['end']}\n")
            f.write(f"{amp['sequence']}\n")

def generate_html_report(set_name, amplicons, primers, extension_time, 
                         polymerase, output_file):
    """Generate HTML validation report"""
    
    # Calculate statistics
    total_amplicons = len(amplicons)
    
    if total_amplicons == 0:
        print(f"Warning: No amplicons generated for {set_name}", file=sys.stderr)
        
        # Write minimal report
        with open(output_file, 'w') as f:
            f.write(f"""<!DOCTYPE html>
<html>
<head>
    <title>Validation Report: {set_name}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .error {{ color: red; font-weight: bold; }}
    </style>
</head>
<body>
    <h1>Validation Report: {set_name}</h1>
    <p class="error">⚠ No amplicons were generated for this primer set!</p>
    <p>This could indicate:</p>
    <ul>
        <li>Primers do not bind to the reference genome</li>
        <li>No valid amplicons within size constraints</li>
        <li>Extension time too short for expected amplicon sizes</li>
    </ul>
</body>
</html>""")
        return
    
    avg_size = sum(a['size'] for a in amplicons) / total_amplicons
    min_size = min(a['size'] for a in amplicons)
    max_size = max(a['size'] for a in amplicons)
    
    # Count amplicons per pair
    pair_counts = defaultdict(int)
    for amp in amplicons:
        pair_counts[amp['pair_id']] += 1
    
    # Count amplicons per contig
    contig_counts = defaultdict(int)
    for amp in amplicons:
        contig_counts[amp['contig']] += 1
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Validation Report: {set_name}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #2c3e50;
        }}
        h2 {{
            color: #34495e;
        }}
        .summary {{
            background-color: white;
            padding: 20px;
            border-radius: 5px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .stat-box {{
            display: inline-block;
            background-color: #3498db;
            color: white;
            padding: 15px;
            margin: 10px;
            border-radius: 5px;
            min-width: 150px;
            text-align: center;
        }}
        .stat-value {{
            font-size: 28px;
            font-weight: bold;
        }}
        .stat-label {{
            font-size: 14px;
            margin-top: 5px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background-color: white;
            margin-top: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: bold;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .primer-seq {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            color: #2c3e50;
        }}
    </style>
</head>
<body>
    <h1>🔬 In Silico PCR Validation Report</h1>
    <h2>{set_name}</h2>
    
    <div class="summary">
        <h3>PCR Conditions</h3>
        <p><strong>Extension Time:</strong> {extension_time} seconds</p>
        <p><strong>Polymerase:</strong> {polymerase}</p>
        <p><strong>Primer Pairs:</strong> {len(primers)}</p>
    </div>
    
    <div class="summary">
        <h3>Amplicon Statistics</h3>
        <div class="stat-box">
            <div class="stat-value">{total_amplicons:,}</div>
            <div class="stat-label">Total Amplicons</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{avg_size:.0f} bp</div>
            <div class="stat-label">Average Size</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{min_size}-{max_size} bp</div>
            <div class="stat-label">Size Range</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{len(contig_counts)}</div>
            <div class="stat-label">Contigs Covered</div>
        </div>
    </div>
    
    <div class="summary">
        <h3>Amplicons per Primer Pair</h3>
        <table>
            <thead>
                <tr>
                    <th>Pair ID</th>
                    <th>Forward Primer (5'→3')</th>
                    <th>Reverse Primer (5'→3')</th>
                    <th>Amplicons</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for pair_id in sorted(pair_counts.keys()):
        primer = primers[pair_id]
        count = pair_counts[pair_id]
        
        html += f"""
                <tr>
                    <td><strong>Pair {pair_id + 1}</strong></td>
                    <td class="primer-seq">{primer['fwd']}</td>
                    <td class="primer-seq">{primer['rev']}</td>
                    <td>{count:,}</td>
                </tr>
"""
    
    html += """
            </tbody>
        </table>
    </div>
    
    <div class="summary">
        <h3>Amplicon Distribution by Contig</h3>
        <table>
            <thead>
                <tr>
                    <th>Contig</th>
                    <th>Amplicons</th>
                    <th>Percentage</th>
                </tr>
            </thead>
            <tbody>
"""
    
    for contig in sorted(contig_counts.keys(), key=lambda x: contig_counts[x], reverse=True):
        count = contig_counts[contig]
        percent = (count / total_amplicons) * 100
        
        html += f"""
                <tr>
                    <td>{contig}</td>
                    <td>{count:,}</td>
                    <td>{percent:.2f}%</td>
                </tr>
"""
    
    html += """
            </tbody>
        </table>
    </div>
    
    <div class="summary">
        <h3>Summary</h3>
        <p>✓ Successfully simulated {total_amplicons:,} amplicons from {len(primers)} primer pairs</p>
        <p>✓ Amplicons span {len(contig_counts)} contig(s)</p>
        <p>✓ Average amplicon size: {avg_size:.0f} bp</p>
    </div>
    
</body>
</html>
""".format(
        total_amplicons=total_amplicons,
        len_primers=len(primers),
        len_contigs=len(contig_counts),
        avg_size=avg_size
    )
    
    with open(output_file, 'w') as f:
        f.write(html)

def generate_summary_report(all_results, output_file):
    """Generate summary report comparing all validated sets"""
    
    html = """<!DOCTYPE html>
<html>
<head>
    <title>Validation Summary Report</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #2c3e50;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            background-color: white;
            margin-top: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #3498db;
            color: white;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
        .rank-1 {
            background-color: #ffd700 !important;
            font-weight: bold;
        }
        .rank-2 {
            background-color: #c0c0c0 !important;
        }
        .rank-3 {
            background-color: #cd7f32 !important;
        }
    </style>
</head>
<body>
    <h1>🧬 Validation Summary Report</h1>
    <p>Comparison of validated primer sets</p>
    
    <table>
        <thead>
            <tr>
                <th>Rank</th>
                <th>Set ID</th>
                <th>Primer Pairs</th>
                <th>Total Amplicons</th>
                <th>Avg Size (bp)</th>
                <th>Contigs</th>
            </tr>
        </thead>
        <tbody>
"""
    
    for rank, result in enumerate(all_results, 1):
        row_class = ''
        if rank == 1:
            row_class = 'rank-1'
        elif rank == 2:
            row_class = 'rank-2'
        elif rank == 3:
            row_class = 'rank-3'
        
        html += f"""
            <tr class="{row_class}">
                <td>{rank}</td>
                <td>{result['set_name']}</td>
                <td>{result['num_pairs']}</td>
                <td>{result['total_amplicons']:,}</td>
                <td>{result['avg_size']:.0f}</td>
                <td>{result['num_contigs']}</td>
            </tr>
"""
    
    html += """
        </tbody>
    </table>
    
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)

def validate_sets(ranked_file, optimized_dir, genome_file, extension_time, 
                 polymerase, output_dir, top_n):
    """
    Validate top N primer sets
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load genome
    print("Loading genome...", file=sys.stderr)
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))
    print(f"✓ Loaded {len(genome)} contigs", file=sys.stderr)
    
    # Get top N sets from ranked file
    print(f"\nParsing ranked sets (top {top_n})...", file=sys.stderr)
    set_ids = parse_ranked_sets(ranked_file, top_n)
    print(f"✓ Found {len(set_ids)} sets to validate", file=sys.stderr)
    
    # Validate each set
    all_results = []
    
    for i, set_id in enumerate(set_ids, 1):
        print(f"\n[{i}/{len(set_ids)}] Validating {set_id}...", file=sys.stderr)
        
        # Find the corresponding primer file
        set_file = os.path.join(optimized_dir, f"{set_id}.primers.txt")
        
        if not os.path.exists(set_file):
            print(f"  ⚠ Warning: File not found: {set_file}", file=sys.stderr)
            continue
        
        # Validate
        amplicons, primers = validate_primer_set(
            set_file, genome, extension_time, polymerase, output_dir
        )
        
        print(f"  ✓ Simulated {len(amplicons)} amplicons from {len(primers)} pairs", 
              file=sys.stderr)
        
        # Generate outputs
        base_name = f"{set_id}_validation"
        
        # BED file
        bed_file = os.path.join(output_dir, f"{base_name}.bed")
        generate_bed_file(amplicons, bed_file)
        
        # FASTA file
        fasta_file = os.path.join(output_dir, f"{base_name}.fasta")
        generate_fasta_file(amplicons, fasta_file)
        
        # HTML report
        html_file = os.path.join(output_dir, f"{base_name}.html")
        generate_html_report(set_id, amplicons, primers, extension_time, 
                           polymerase, html_file)
        
        # Collect results for summary
        if amplicons:
            contig_counts = defaultdict(int)
            for amp in amplicons:
                contig_counts[amp['contig']] += 1
            
            avg_size = sum(a['size'] for a in amplicons) / len(amplicons)
            
            all_results.append({
                'set_name': set_id,
                'num_pairs': len(primers),
                'total_amplicons': len(amplicons),
                'avg_size': avg_size,
                'num_contigs': len(contig_counts)
            })
    
    # Generate summary report
    if all_results:
        summary_file = os.path.join(output_dir, "validation_summary.html")
        generate_summary_report(all_results, summary_file)
        print(f"\n✓ Summary report: {summary_file}", file=sys.stderr)
    
    print(f"\n✓ Validation complete! Results in: {output_dir}/", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Validate primer sets through in silico PCR simulation'
    )
    parser.add_argument('--ranked-sets', required=True, 
                       help='Ranked sets TSV file')
    parser.add_argument('--optimized-dir', required=True,
                       help='Directory with optimized primer set files')
    parser.add_argument('--genome', required=True, help='Genome FASTA')
    parser.add_argument('--extension-time', type=int, required=True,
                       help='PCR extension time (seconds)')
    parser.add_argument('--polymerase', default='Q5', 
                       help='Polymerase (Q5 or Taq)')
    parser.add_argument('--output-dir', required=True, 
                       help='Output directory for validation results')
    parser.add_argument('--top-n', type=int, default=10,
                       help='Validate top N sets')
    
    args = parser.parse_args()
    
    validate_sets(
        ranked_file=args.ranked_sets,
        optimized_dir=args.optimized_dir,
        genome_file=args.genome,
        extension_time=args.extension_time,
        polymerase=args.polymerase,
        output_dir=args.output_dir,
        top_n=args.top_n
    )

if __name__ == '__main__':
    main()
