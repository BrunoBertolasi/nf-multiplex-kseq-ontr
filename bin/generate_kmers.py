#!/usr/bin/env python3
"""
Generate candidate primers from target sequences
Uses k-mer sliding window approach
"""

import sys
import argparse
try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False
    print("Warning: primer3-py not installed, using basic Tm calculation", file=sys.stderr)

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import MeltingTemp as mt

def calc_tm(seq):
    """Calculate melting temperature"""
    if HAS_PRIMER3:
        return primer3.calc_tm(seq)
    else:
        # Basic Tm calculation
        return mt.Tm_NN(seq)

def is_valid_primer(seq, min_gc=40, max_gc=60, min_tm=55, max_tm=65, max_homopolymer=4):
    """Check if primer meets basic criteria"""
    
    # Convert to string and uppercase
    seq = str(seq).upper()
    
    # No ambiguous bases
    if any(b not in 'ACGT' for b in seq):
        return False
    
    # GC content
    gc = gc_fraction(seq) * 100
    if not (min_gc <= gc <= max_gc):
        return False
    
    # Homopolymer check
    for nt in 'ACGT':
        if nt * (max_homopolymer + 1) in seq:
            return False
    
    # Tm calculation
    try:
        tm = calc_tm(seq)
        if not (min_tm <= tm <= max_tm):
            return False
    except:
        # If Tm calculation fails, skip this primer
        return False
    
    return True

def generate_candidates(target_fasta, primer_length, output_forward, output_reverse, **kwargs):
    """
    Generate all valid k-mer primers from target sequences
    """
    
    forward_primers = {}
    reverse_primers = {}
    
    for record in SeqIO.parse(target_fasta, 'fasta'):
        seq = str(record.seq).upper()
        
        # Sliding window for forward primers
        for i in range(len(seq) - primer_length + 1):
            kmer = seq[i:i + primer_length]
            
            if is_valid_primer(kmer, **kwargs):
                if kmer not in forward_primers:
                    forward_primers[kmer] = {
                        'tm': calc_tm(kmer),
                        'gc': gc_fraction(kmer) * 100
                    }
        
        # Reverse complement for reverse primers
        rc_seq = str(record.seq.reverse_complement()).upper()
        for i in range(len(rc_seq) - primer_length + 1):
            kmer = rc_seq[i:i + primer_length]
            
            if is_valid_primer(kmer, **kwargs):
                if kmer not in reverse_primers:
                    reverse_primers[kmer] = {
                        'tm': calc_tm(kmer),
                        'gc': gc_fraction(kmer) * 100
                    }
    
    # Write output
    with open(output_forward, 'w') as f:
        for primer in sorted(forward_primers.keys()):
            info = forward_primers[primer]
            f.write(f"{primer}\t{info['tm']:.1f}\t{info['gc']:.1f}\n")
    
    with open(output_reverse, 'w') as f:
        for primer in sorted(reverse_primers.keys()):
            info = reverse_primers[primer]
            f.write(f"{primer}\t{info['tm']:.1f}\t{info['gc']:.1f}\n")
    
    print(f"Generated {len(forward_primers)} forward primers", file=sys.stderr)
    print(f"Generated {len(reverse_primers)} reverse primers", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Generate candidate primers from target sequences'
    )
    parser.add_argument('--input', required=True, help='Input FASTA file')
    parser.add_argument('--primer-length', type=int, default=18, help='Primer length')
    parser.add_argument('--min-gc', type=float, default=40, help='Minimum GC%')
    parser.add_argument('--max-gc', type=float, default=60, help='Maximum GC%')
    parser.add_argument('--min-tm', type=float, default=55, help='Minimum Tm (°C)')
    parser.add_argument('--max-tm', type=float, default=65, help='Maximum Tm (°C)')
    parser.add_argument('--max-homopolymer', type=int, default=4, help='Max homopolymer length')
    parser.add_argument('--output-forward', required=True, help='Output forward primers')
    parser.add_argument('--output-reverse', required=True, help='Output reverse primers')
    
    args = parser.parse_args()
    
    generate_candidates(
        target_fasta=args.input,
        primer_length=args.primer_length,
        output_forward=args.output_forward,
        output_reverse=args.output_reverse,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        min_tm=args.min_tm,
        max_tm=args.max_tm,
        max_homopolymer=args.max_homopolymer
    )

if __name__ == '__main__':
    main()
