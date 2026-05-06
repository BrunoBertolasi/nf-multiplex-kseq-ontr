#!/usr/bin/env python3
"""
Additional utility for checking primer dimers and hairpins
Can be used as a standalone validation tool
"""

import sys
import argparse
from Bio.Seq import Seq

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False
    print("Warning: primer3-py not installed, using basic checks", file=sys.stderr)

def reverse_complement(seq):
    """Get reverse complement of sequence"""
    return str(Seq(seq).reverse_complement())

def check_self_dimer(seq, max_3prime_matches=2):
    """Check if primer forms self-dimers"""
    rc = reverse_complement(seq)
    
    # Check 3' end (last 5 bases)
    end = seq[-5:]
    rc_end = rc[-5:]
    
    matches = sum(a == b for a, b in zip(end, rc_end))
    
    if matches > max_3prime_matches:
        return True, f"Self-dimer: {matches} 3' matches"
    
    return False, "OK"

def check_hetero_dimer(seq1, seq2, max_3prime_matches=2, max_internal_matches=6):
    """Check if two primers form heterodimers"""
    rc2 = reverse_complement(seq2)
    
    # Check 3' complementarity
    end1 = seq1[-5:]
    end2_rc = rc2[-5:]
    
    end_matches = sum(a == b for a, b in zip(end1, end2_rc))
    
    if end_matches > max_3prime_matches:
        return True, f"Heterodimer: {end_matches} 3' matches"
    
    # Check internal complementarity
    max_continuous = 0
    for i in range(len(seq1) - 4):
        for j in range(len(rc2) - 4):
            window_size = min(5, len(seq1) - i, len(rc2) - j)
            matches = sum(seq1[i + k] == rc2[j + k] for k in range(window_size))
            max_continuous = max(max_continuous, matches)
    
    if max_continuous > max_internal_matches:
        return True, f"Internal complementarity: {max_continuous} continuous matches"
    
    return False, "OK"

def check_hairpin(seq, max_hairpin_tm=45):
    """Check for hairpin formation"""
    if HAS_PRIMER3:
        try:
            tm = primer3.calc_hairpin_tm(seq)
            if tm > max_hairpin_tm:
                return True, f"Hairpin Tm: {tm:.1f}°C"
        except:
            pass
    
    # Basic hairpin check (palindromes)
    rc = reverse_complement(seq)
    
    for i in range(len(seq) - 5):
        for j in range(len(rc) - 5):
            # Check for 6+ bp complementarity within same primer
            if seq[i:i+6] == rc[j:j+6]:
                return True, "Potential hairpin: 6bp palindrome detected"
    
    return False, "OK"

def analyze_primer_set(primer_file, output_file, report_file):
    """
    Analyze a primer set for potential issues
    """
    
    # Read primers
    primers = []
    with open(primer_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                primers.append({
                    'fwd': parts[0],
                    'rev': parts[1]
                })
    
    print(f"Analyzing {len(primers)} primer pairs...", file=sys.stderr)
    
    # Check for issues
    issues = []
    
    for i, pair in enumerate(primers):
        fwd = pair['fwd']
        rev = pair['rev']
        
        # Self-dimer checks
        fwd_self, fwd_msg = check_self_dimer(fwd)
        if fwd_self:
            issues.append(f"Pair {i+1} FWD: {fwd_msg}")
        
        rev_self, rev_msg = check_self_dimer(rev)
        if rev_self:
            issues.append(f"Pair {i+1} REV: {rev_msg}")
        
        # Heterodimer check
        hetero, hetero_msg = check_hetero_dimer(fwd, rev)
        if hetero:
            issues.append(f"Pair {i+1} FWD-REV: {hetero_msg}")
        
        # Hairpin checks
        fwd_hp, fwd_hp_msg = check_hairpin(fwd)
        if fwd_hp:
            issues.append(f"Pair {i+1} FWD: {fwd_hp_msg}")
        
        rev_hp, rev_hp_msg = check_hairpin(rev)
        if rev_hp:
            issues.append(f"Pair {i+1} REV: {rev_hp_msg}")
        
        # Check against all other primers in set
        for j, other_pair in enumerate(primers):
            if i >= j:
                continue
            
            # Check all combinations
            for primer_a, name_a in [(fwd, 'FWD'), (rev, 'REV')]:
                for primer_b, name_b in [(other_pair['fwd'], 'FWD'), (other_pair['rev'], 'REV')]:
                    hetero, msg = check_hetero_dimer(primer_a, primer_b)
                    if hetero:
                        issues.append(f"Pair {i+1}-{name_a} × Pair {j+1}-{name_b}: {msg}")
    
    # Write issues to file
    with open(output_file, 'w') as f:
        if issues:
            f.write(f"Found {len(issues)} potential issues:\n\n")
            for issue in issues:
                f.write(f"{issue}\n")
        else:
            f.write("✓ No dimer or hairpin issues detected!\n")
    
    # Generate report
    with open(report_file, 'w') as f:
        f.write("Primer Dimer & Hairpin Analysis Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total primer pairs analyzed: {len(primers)}\n")
        f.write(f"Potential issues found: {len(issues)}\n\n")
        
        if issues:
            f.write("Issues:\n")
            f.write("-" * 50 + "\n")
            for issue in issues:
                f.write(f"  • {issue}\n")
        else:
            f.write("✓ All primers passed dimer and hairpin checks!\n")
    
    print(f"\n{'⚠' if issues else '✓'} Analysis complete", file=sys.stderr)
    print(f"  Issues found: {len(issues)}", file=sys.stderr)
    print(f"  Results: {output_file}", file=sys.stderr)
    print(f"  Report: {report_file}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Check primers for dimers and hairpins'
    )
    parser.add_argument('--primers', required=True, 
                       help='Primer file (TSV with FWD and REV columns)')
    parser.add_argument('--output', required=True, help='Output issues file')
    parser.add_argument('--report', required=True, help='Output report file')
    parser.add_argument('--max-3prime', type=int, default=2,
                       help='Maximum 3\' complementarity')
    parser.add_argument('--max-internal', type=int, default=6,
                       help='Maximum internal complementarity')
    parser.add_argument('--max-hairpin-tm', type=float, default=45,
                       help='Maximum hairpin Tm (°C)')
    
    args = parser.parse_args()
    
    analyze_primer_set(args.primers, args.output, args.report)

if __name__ == '__main__':
    main()
