#!/usr/bin/env python3
"""
Count primer binding sites across genome with fuzzy matching
Applies filters for binding site count and target region specificity
"""

import sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

# Try to import regex for fuzzy matching
try:
    import regex
    HAS_REGEX = True
except ImportError:
    HAS_REGEX = False
    print("Warning: regex module not installed, using exact matching only", file=sys.stderr)

def fuzzy_search(pattern, text, max_mismatches=1):
    """
    Find all occurrences of pattern in text allowing mismatches
    Returns list of (position, matched_sequence) tuples
    """
    if HAS_REGEX and max_mismatches > 0:
        # Use regex for fuzzy matching
        fuzzy_pattern = f"({pattern}){{s<={max_mismatches}}}"
        matches = []
        
        for match in regex.finditer(fuzzy_pattern, text, overlapped=True):
            matches.append((match.start(), match.group()))
        
        return matches
    else:
        # Exact matching only
        matches = []
        pos = 0
        while True:
            pos = text.find(pattern, pos)
            if pos == -1:
                break
            matches.append((pos, pattern))
            pos += 1
        return matches

def load_target_regions(bed_file):
    """
    Load target regions from BED file
    Returns dict: {contig: [(start, end), ...]}
    """
    regions = defaultdict(list)
    
    try:
        with open(bed_file) as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 3:
                        contig = fields[0]
                        start = int(fields[1])
                        end = int(fields[2])
                        regions[contig].append((start, end))
    except FileNotFoundError:
        pass  # No target regions file
    
    return regions

def is_in_target_region(contig, pos, target_regions):
    """Check if position is within any target region"""
    if not target_regions:
        return True  # No restrictions
    
    if contig not in target_regions:
        return False
    
    for start, end in target_regions[contig]:
        if start <= pos <= end:
            return True
    
    return False

def count_binding_sites(primer_file, genome_file, target_bed, max_mismatches, 
                       strict_mode, min_count, max_count, output):
    """
    Count how many times each primer binds to genome
    Filter by target regions if strict_mode enabled
    """
    
    # Load genome
    print("Loading genome...", file=sys.stderr)
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))
    print(f"Loaded {len(genome)} contigs", file=sys.stderr)
    
    # Load target regions (if strict mode)
    target_regions = load_target_regions(target_bed) if strict_mode else {}
    
    if strict_mode and target_regions:
        total_regions = sum(len(v) for v in target_regions.values())
        print(f"Strict mode: restricting to {total_regions} target regions", file=sys.stderr)
    
    # Process primers
    total_primers = 0
    valid_primers = 0
    
    with open(primer_file) as f:
        primers = [line.strip().split('\t') for line in f if line.strip()]
    
    print(f"Processing {len(primers)} primers...", file=sys.stderr)
    
    with open(output, 'w') as out:
        for idx, primer_data in enumerate(primers):
            if len(primer_data) < 3:
                continue
            
            primer, tm, gc = primer_data[:3]
            total_primers += 1
            
            if (idx + 1) % 1000 == 0:
                print(f"  Processed {idx + 1}/{len(primers)} primers...", file=sys.stderr)
            
            # Find all binding sites
            sites = []
            
            for contig_id, record in genome.items():
                seq = str(record.seq).upper()
                
                # Search forward strand
                fwd_matches = fuzzy_search(primer, seq, max_mismatches)
                for pos, matched in fwd_matches:
                    # Check if in target region (if strict mode)
                    if not strict_mode or is_in_target_region(contig_id, pos, target_regions):
                        sites.append({
                            'contig': contig_id,
                            'pos': pos,
                            'strand': '+',
                            'match': matched
                        })
                
                # Search reverse strand (reverse complement)
                primer_rc = str(Seq(primer).reverse_complement())
                rev_matches = fuzzy_search(primer_rc, seq, max_mismatches)
                for pos, matched in rev_matches:
                    if not strict_mode or is_in_target_region(contig_id, pos, target_regions):
                        sites.append({
                            'contig': contig_id,
                            'pos': pos,
                            'strand': '-',
                            'match': matched
                        })
            
            # Filter by site count
            site_count = len(sites)
            
            if min_count <= site_count <= max_count:
                valid_primers += 1
                
                # Write primer with binding site count
                out.write(f"{primer}\t{tm}\t{gc}\t{site_count}\t")
                
                # Write first 20 sites (for reference)
                site_strings = [
                    f"{s['contig']}:{s['pos']}:{s['strand']}" 
                    for s in sites[:20]
                ]
                out.write(",".join(site_strings))
                out.write("\n")
    
    print(f"\nTotal primers: {total_primers}", file=sys.stderr)
    print(f"Valid primers (count {min_count}-{max_count}): {valid_primers}", file=sys.stderr)
    print(f"Filtered out: {total_primers - valid_primers}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Count primer binding sites across genome'
    )
    parser.add_argument('--primers', required=True, help='Primer file')
    parser.add_argument('--genome', required=True, help='Genome FASTA')
    parser.add_argument('--target-bed', required=True, help='Target regions BED file')
    parser.add_argument('--max-mismatches', type=int, default=1, 
                       help='Maximum mismatches allowed')
    parser.add_argument('--strict-mode', action='store_true',
                       help='Enforce target region restrictions')
    parser.add_argument('--min-count', type=int, default=100,
                       help='Minimum binding sites')
    parser.add_argument('--max-count', type=int, default=1000,
                       help='Maximum binding sites')
    parser.add_argument('--output', required=True, help='Output file')
    
    args = parser.parse_args()
    
    count_binding_sites(
        primer_file=args.primers,
        genome_file=args.genome,
        target_bed=args.target_bed,
        max_mismatches=args.max_mismatches,
        strict_mode=args.strict_mode,
        min_count=args.min_count,
        max_count=args.max_count,
        output=args.output
    )

if __name__ == '__main__':
    main()
