#!/usr/bin/env python3
"""
Find maximum compatible sets of primer pairs
Uses graph-based approach to identify non-interacting primers
"""

import sys
import argparse
import os
from Bio.Seq import Seq

# Try to import libraries
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    print("Warning: networkx not installed, using simple algorithm", file=sys.stderr)

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False

def parse_primer_pairs(pair_file):
    """Parse primer pairs file"""
    pairs = []
    
    with open(pair_file) as f:
        for line in f:
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

def check_3prime_complementarity(seq1, seq2, max_matches=2):
    """
    Check 3' end complementarity (primer dimer risk)
    Returns True if dimers likely to form
    """
    # Check last 5 bases
    end1 = seq1[-5:]
    end2_rc = str(Seq(seq2[-5:]).reverse_complement())
    
    matches = sum(a == b for a, b in zip(end1, end2_rc))
    
    return matches > max_matches

def check_internal_complementarity(seq1, seq2, max_matches=6):
    """
    Check internal complementarity (cross-priming risk)
    """
    seq2_rc = str(Seq(seq2).reverse_complement())
    
    # Count maximum continuous matches in sliding window
    max_continuous = 0
    
    for i in range(len(seq1) - 4):
        for j in range(len(seq2_rc) - 4):
            window_size = min(5, len(seq1) - i, len(seq2_rc) - j)
            matches = sum(
                seq1[i + k] == seq2_rc[j + k] 
                for k in range(window_size)
            )
            max_continuous = max(max_continuous, matches)
    
    return max_continuous > max_matches

def check_hairpin(seq, max_tm=45):
    """Check for hairpin formation"""
    if HAS_PRIMER3:
        try:
            result = primer3.calc_hairpin_tm(seq)
            return result > max_tm
        except:
            return False
    return False

def primers_interact(pair1, pair2, max_3prime=2, max_internal=6):
    """
    Check if two primer pairs interact
    Returns True if incompatible
    """
    # All primers from both pairs
    primers = [
        pair1['fwd'],
        pair1['rev'],
        pair2['fwd'],
        pair2['rev']
    ]
    
    # Check all combinations
    for i, p1 in enumerate(primers):
        for p2 in primers[i + 1:]:
            # 3' complementarity check
            if check_3prime_complementarity(p1, p2, max_3prime):
                return True
            
            # Internal complementarity check
            if check_internal_complementarity(p1, p2, max_internal):
                return True
    
    return False

def find_compatible_sets_networkx(pairs, min_pairs, max_pairs, max_3prime, max_internal):
    """Find compatible sets using NetworkX"""
    
    # Build incompatibility graph
    G = nx.Graph()
    
    for i in range(len(pairs)):
        G.add_node(i)
    
    print("Checking primer interactions...", file=sys.stderr)
    interaction_count = 0
    
    for i in range(len(pairs)):
        if (i + 1) % 10 == 0:
            print(f"  Checked {i + 1}/{len(pairs)} pairs", file=sys.stderr)
        
        for j in range(i + 1, len(pairs)):
            if primers_interact(pairs[i], pairs[j], max_3prime, max_internal):
                G.add_edge(i, j)
                interaction_count += 1
    
    print(f"Found {interaction_count} incompatible primer pair interactions", file=sys.stderr)
    
    # Find maximum independent sets (compatible primers)
    complement = nx.complement(G)
    
    print("Finding maximum compatible sets...", file=sys.stderr)
    cliques = list(nx.find_cliques(complement))
    
    # Sort by size (largest first)
    cliques.sort(key=len, reverse=True)
    
    # Filter by size requirements
    valid_sets = [c for c in cliques if min_pairs <= len(c) <= max_pairs]
    
    return valid_sets, G

def find_compatible_sets_simple(pairs, min_pairs, max_pairs, max_3prime, max_internal):
    """Simple greedy algorithm without NetworkX"""
    
    print("Using simple greedy algorithm (NetworkX not available)", file=sys.stderr)
    
    # Build incompatibility matrix
    incompatible = [[False] * len(pairs) for _ in range(len(pairs))]
    
    for i in range(len(pairs)):
        if (i + 1) % 10 == 0:
            print(f"  Checked {i + 1}/{len(pairs)} pairs", file=sys.stderr)
        
        for j in range(i + 1, len(pairs)):
            if primers_interact(pairs[i], pairs[j], max_3prime, max_internal):
                incompatible[i][j] = True
                incompatible[j][i] = True
    
    # Greedy set building
    sets = []
    remaining = set(range(len(pairs)))
    
    while remaining and len(sets) < 100:  # Limit to 100 sets
        current_set = []
        available = remaining.copy()
        
        while available:
            # Pick primer with fewest incompatibilities
            best = min(available, key=lambda x: sum(incompatible[x][y] for y in available))
            current_set.append(best)
            
            # Remove incompatible primers
            available = {p for p in available if not incompatible[best][p] and p != best}
        
        if min_pairs <= len(current_set) <= max_pairs:
            sets.append(current_set)
        
        # Remove used primers for next iteration
        if current_set:
            remaining.discard(current_set[0])
    
    return sets, None

def find_compatible_sets(pair_file, min_pairs, max_pairs, max_3prime, 
                        max_internal, output_dir, graph_file):
    """
    Find all maximal compatible primer sets
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse pairs
    print("Loading primer pairs...", file=sys.stderr)
    pairs = parse_primer_pairs(pair_file)
    print(f"Total primer pairs: {len(pairs)}", file=sys.stderr)
    
    # Find compatible sets
    if HAS_NETWORKX:
        valid_sets, graph = find_compatible_sets_networkx(
            pairs, min_pairs, max_pairs, max_3prime, max_internal
        )
    else:
        valid_sets, graph = find_compatible_sets_simple(
            pairs, min_pairs, max_pairs, max_3prime, max_internal
        )
    
    print(f"\nFound {len(valid_sets)} compatible sets", file=sys.stderr)
    
    # Write graph info
    if graph and HAS_NETWORKX:
        with open(graph_file, 'w') as f:
            f.write(f"Nodes (primer pairs): {graph.number_of_nodes()}\n")
            f.write(f"Edges (incompatibilities): {graph.number_of_edges()}\n")
            f.write(f"Compatible sets found: {len(valid_sets)}\n")
    
    # Write sets (all of them, sorted by size)
    valid_sets.sort(key=len, reverse=True)
    
    for set_id, primer_set in enumerate(valid_sets):
        with open(f"{output_dir}/set_{set_id:04d}.txt", 'w') as f:
            f.write(f"# Compatible primer set {set_id}\n")
            f.write(f"# Size: {len(primer_set)} pairs\n")
            f.write("# Forward\tReverse\tFwd_Tm\tRev_Tm\tAmplicons\tAvg_Size\tMin_Size\tMax_Size\n")
            
            for idx in primer_set:
                pair = pairs[idx]
                f.write(f"{pair['fwd']}\t{pair['rev']}\t")
                f.write(f"{pair['fwd_tm']:.1f}\t{pair['rev_tm']:.1f}\t")
                f.write(f"{pair['count']}\t{pair['avg_size']:.0f}\t")
                f.write(f"{pair['min_size']}\t{pair['max_size']}\n")
    
    print(f"\n✓ Written {len(valid_sets)} compatible sets to {output_dir}/", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Find compatible primer sets'
    )
    parser.add_argument('--primer-pairs', required=True, help='Primer pairs file')
    parser.add_argument('--min-pairs', type=int, default=15, help='Minimum pairs per set')
    parser.add_argument('--max-pairs', type=int, default=50, help='Maximum pairs per set')
    parser.add_argument('--max-3prime', type=int, default=2, 
                       help='Max 3\' complementarity')
    parser.add_argument('--max-internal', type=int, default=6,
                       help='Max internal complementarity')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--graph', required=True, help='Graph info output file')
    
    args = parser.parse_args()
    
    find_compatible_sets(
        pair_file=args.primer_pairs,
        min_pairs=args.min_pairs,
        max_pairs=args.max_pairs,
        max_3prime=args.max_3prime,
        max_internal=args.max_internal,
        output_dir=args.output_dir,
        graph_file=args.graph
    )

if __name__ == '__main__':
    main()
