# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    species = {
        "Gallus gallus": gg_seq,
        "Mus musculus": mm_seq,
        "Balaeniceps rex": br_seq,
        "Tursiops truncatus": tt_seq
    }
    
    results = {}

    # Align all species to humans and print species in order of most similar to human BRD
    for name, seq in species.items():
        nw = NeedlemanWunsch(
            'substitution_matrices/BLOSUM62.mat',
            gap_open=-10,
            gap_extend=-1
        )
        score,seq_align_ref, seq_align = nw.align(hs_seq, seq)
        results[name] = {
            'ref':seq_align_ref,
            'aln':seq_align,
            'score':score
        }
    
    sorted_species = dict(sorted(results.items(), key=lambda x: x[1]['score'], reverse=True))

    print("Species in order of most similar to human BRD.")
    
    print("\n\tAlignment scores:\n")
    for name, data in sorted_species.items():
        print(f"\t{name}: {data['score']}")


if __name__ == "__main__":
    main()
