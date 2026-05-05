#!/usr/bin/env python3
"""
Pair forward and reverse primers to find valid amplicons
Applies extension time constraint to filter amplicon sizes
"""

import sys
import argparse
from collections import defaultdict
from Bio import SeqIO

def calculate_max_amplicon_size(extension_time, polymerase='Q5'):
    """
    Calculate maximum amplicon size based on extension time
    
    Q5: 15-30 sec/kb (we use conservative 30 sec/kb)
    Taq: 60 sec/kb
    """
    if polymerase.upper() == 'Q5':
        # Conservative: 30 sec/kb
        max_size = int(extension_time / 30 * 1000)
    else:  # Taq
        max_size = int(extension_time / 60 * 1000)
    
    return max_size

def parse_primer_sites(primer_file):
    """
    Parse primer file with binding sites
    Format: primer\ttm\tgc\tcount\tsites
    """
    primers = []
    
    with open(primer_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            primer, tm, gc, count, sites_str = parts[:5]
            
            # Parse sites
            sites = []
            if sites_str:
                for site in sites_str.split(','):
                    try:
                        contig, pos, strand = site.split(':')
                        sites.append({
                            'contig': contig,
                            'pos': int(pos),
                            'strand': strand
                        })
                    except ValueError:
                        continue
            
            primers.append({
                'seq': primer,
                'tm': float(tm),
                'gc': float(gc),
                'count': int(count),
                'sites': sites
            })
    
    return primers

def pair_primers(forward_file, reverse_file, genome, genome_index, 
                min_size, max_size, extension_time, polymerase,
                min_amplicons, max_amplicons, output, stats_file):
    """
    Find valid F/R primer pairs considering extension time constraint
    """
    
    # Calculate actual max size from extension time
    extension_max = calculate_max_amplicon_size(extension_time, polymerase)
    actual_max = min(max_size, extension_max)
    
    print(f"Extension time: {extension_time}s with {polymerase}", file=sys.stderr)
    print(f"Maximum amplicon size from extension: {extension_max} bp", file=sys.stderr)
    print(f"Actual maximum amplicon size: {actual_max} bp", file=sys.stderr)
    
    # Parse primers
    print("\nLoading primers...", file=sys.stderr)
    forward_primers = parse_primer_sites(forward_file)
    reverse_primers = parse_primer_sites(reverse_file)
    
    print(f"Forward primers: {len(forward_primers)}", file=sys.stderr)
    print(f"Reverse primers: {len(reverse_primers)}", file=sys.stderr)
    
    # Index forward sites by contig for fast lookup
    print("\nIndexing forward primer sites...", file=sys.stderr)
    fwd_by_contig = defaultdict(list)
    for i, fwd in enumerate(forward_primers):
        for site in fwd['sites']:
            fwd_by_contig[site['contig']].append({
                'primer_idx': i,
                'pos': site['pos'],
                'strand': site['strand']
            })
    
    # Find pairs
    print("\nFinding primer pairs...", file=sys.stderr)
    pairs = []
    
    for rev_idx, rev in enumerate(reverse_primers):
        if (rev_idx + 1) % 100 == 0:
            print(f"  Processed {rev_idx + 1}/{len(reverse_primers)} reverse primers...", 
                  file=sys.stderr)
        
        for rev_site in rev['sites']:
            contig = rev_site['contig']
            rev_pos = rev_site['pos']
            
            # Find forward primers on same contig
            if contig not in fwd_by_contig:
                continue
            
            for fwd_info in fwd_by_contig[contig]:
                fwd_idx = fwd_info['primer_idx']
                fwd_pos = fwd_info['pos']
                
                # Calculate amplicon size
                amplicon_size = abs(rev_pos - fwd_pos)
                
                # Check if within size range AND extension time limit
                if min_size <= amplicon_size <= actual_max:
                    pairs.append({
                        'fwd_idx': fwd_idx,
                        'rev_idx': rev_idx,
                        'contig': contig,
                        'fwd_pos': fwd_pos,
                        'rev_pos': rev_pos,
                        'size': amplicon_size
                    })
    
    print(f"\nFound {len(pairs)} individual amplicons", file=sys.stderr)
    
    # Group by primer pair
    print("Grouping by primer pair...", file=sys.stderr)
    pair_groups = defaultdict(list)
    for pair in pairs:
        key = (pair['fwd_idx'], pair['rev_idx'])
        pair_groups[key].append(pair)
    
    print(f"Total primer pairs: {len(pair_groups)}", file=sys.stderr)
    
    # Filter by amplicon count and write output
    print(f"\nFiltering pairs (amplicon count {min_amplicons}-{max_amplicons})...", 
          file=sys.stderr)
    
    valid_pairs = 0
    
    with open(output, 'w') as out:
        for (fwd_idx, rev_idx), amplicons in pair_groups.items():
            amplicon_count = len(amplicons)
            
            # Filter by amplicon count
            if not (min_amplicons <= amplicon_count <= max_amplicons):
                continue
            
            valid_pairs += 1
            
            fwd = forward_primers[fwd_idx]
            rev = reverse_primers[rev_idx]
            
            avg_size = sum(a['size'] for a in amplicons) / len(amplicons)
            min_amp_size = min(a['size'] for a in amplicons)
            max_amp_size = max(a['size'] for a in amplicons)
            
            out.write(f"{fwd['seq']}\t{rev['seq']}\t")
            out.write(f"{fwd['tm']:.1f}\t{rev['tm']:.1f}\t")
            out.write(f"{amplicon_count}\t{avg_size:.0f}\t{min_amp_size}\t{max_amp_size}\t")
            
            # Write sample amplicon locations (first 20)
            amp_strings = [
                f"{a['contig']}:{a['fwd_pos']}-{a['rev_pos']}" 
                for a in amplicons[:20]
            ]
            out.write(",".join(amp_strings))
            out.write("\n")
    
    # Write statistics
    with open(stats_file, 'w') as stats:
        stats.write("Pairing Statistics\n")
        stats.write("==================\n\n")
        stats.write(f"Forward primers: {len(forward_primers)}\n")
        stats.write(f"Reverse primers: {len(reverse_primers)}\n")
        stats.write(f"Extension time: {extension_time}s ({polymerase})\n")
        stats.write(f"Max amplicon size (extension limit): {extension_max} bp\n")
        stats.write(f"Actual max amplicon size: {actual_max} bp\n")
        stats.write(f"Amplicon size range: {min_size}-{actual_max} bp\n\n")
        stats.write(f"Total amplicons found: {len(pairs)}\n")
        stats.write(f"Total primer pairs: {len(pair_groups)}\n")
        stats.write(f"Valid pairs (count {min_amplicons}-{max_amplicons}): {valid_pairs}\n")
    
    print(f"\n✓ Valid primer pairs: {valid_pairs}", file=sys.stderr)
    print(f"✓ Statistics written to: {stats_file}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Pair forward and reverse primers'
    )
    parser.add_argument('--forward', required=True, help='Forward primers file')
    parser.add_argument('--reverse', required=True, help='Reverse primers file')
    parser.add_argument('--genome', required=True, help='Genome FASTA')
    parser.add_argument('--genome-index', required=True, help='Genome index (.fai)')
    parser.add_argument('--min-size', type=int, required=True, help='Min amplicon size')
    parser.add_argument('--max-size', type=int, required=True, help='Max amplicon size')
    parser.add_argument('--extension-time', type=int, required=True, 
                       help='PCR extension time (seconds)')
    parser.add_argument('--polymerase', default='Q5', help='Polymerase (Q5 or Taq)')
    parser.add_argument('--min-amplicons', type=int, default=100,
                       help='Minimum amplicons per pair')
    parser.add_argument('--max-amplicons', type=int, default=1000,
                       help='Maximum amplicons per pair')
    parser.add_argument('--output', required=True, help='Output primer pairs file')
    parser.add_argument('--stats', required=True, help='Output statistics file')
    
    args = parser.parse_args()
    
    pair_primers(
        forward_file=args.forward,
        reverse_file=args.reverse,
        genome=args.genome,
        genome_index=args.genome_index,
        min_size=args.min_size,
        max_size=args.max_size,
        extension_time=args.extension_time,
        polymerase=args.polymerase,
        min_amplicons=args.min_amplicons,
        max_amplicons=args.max_amplicons,
        output=args.output,
        stats_file=args.stats
    )

if __name__ == '__main__':
    main()
