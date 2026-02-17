import random
import typing

nt_to_aa = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


synonymous_codons = {
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "Q": ["CAA", "CAG"],
    "M": ["ATG"],
    "N": ["AAC", "AAT"],
    "P": ["CCT", "CCG", "CCA", "CCC"],
    "K": ["AAG", "AAA"],
    "_": ["TAG", "TGA", "TAA"],
    "T": ["ACC", "ACA", "ACG", "ACT"],
    "F": ["TTT", "TTC"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "G": ["GGT", "GGG", "GGA", "GGC"],
    "I": ["ATC", "ATA", "ATT"],
    "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "H": ["CAT", "CAC"],
    "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "W": ["TGG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "E": ["GAG", "GAA"],
    "Y": ["TAT", "TAC"],
}

DEFAULT_ORGANISM_PRIORITY = 50


def random_synonymous_codon(codon: str) -> str:
    possible_codons = synonymous_codons[nt_to_aa[codon]]
    return random.choice(possible_codons)


def synonymous_codon_permutation(seq: str, skipped_codons_num: int = 0) -> str:
    if len(seq) % 3 != 0:
        raise ValueError(f"len of seq {seq} is {len(seq)} which is not divisible by 3")
    permutation = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        new_codon = codon
        if i >= skipped_codons_num * 3:
            new_codon = random_synonymous_codon(codon)
        permutation += new_codon
    return permutation


def write_fasta(fid, list_seq, list_name):
    ofile = open(fid + '.fasta', "w+")
    for i in range(len(list_seq)):
        ofile.write(">" + str(list_name[i]) + "\n" + str(list_seq[i]) + "\n")
    ofile.close()


# write ideas for the promoter model
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def translate(seq, table=nt_to_aa):  # DONE
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
        return protein
    else:
        return ValueError('len(seq)%3 !=0')


def change_all_codons_of_aa(
        seq: str,
        selected_codon: str,
        skipped_codons_num: int = 0,
        default_codon: typing.Optional[str] = None,
        should_dedup_codons: bool = False,
) -> str:
    split_seq = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
    new_split_seq = []
    selected_codon_aa = nt_to_aa[selected_codon]
    for i, codon in enumerate(split_seq):
        if i < skipped_codons_num or nt_to_aa[codon] != selected_codon_aa:
            new_split_seq.append(codon)
        elif should_dedup_codons is False:
            new_split_seq.append(selected_codon)
        else:
            should_use_selected_codon = True
            for j in range(i - 1, skipped_codons_num, -1):
                if nt_to_aa[split_seq[j]] != selected_codon_aa:
                    break
                should_use_selected_codon = i - j % 2 == 0
            if should_use_selected_codon:
                new_split_seq.append(selected_codon)
            elif default_codon is not None:
                new_split_seq.append(default_codon)
            else:
                new_split_seq.append(codon)
    return ''.join(new_split_seq)


def unique(list1):
    return sorted(set(list1))


def validate_module_output(original_sequence: str, new_sequence: str) -> None:
    if len(original_sequence) != len(new_sequence):
        raise RuntimeError(f"Optimized sequence length {len(new_sequence)} is different than the initial "
                           f"sequence length {len(original_sequence)}")

    translated_original_seq = translate(original_sequence)
    translated_optimized_seq = translate(new_sequence)
    if translated_original_seq != translated_optimized_seq:
        raise RuntimeError(f"Optimized sequence is translated to {translated_optimized_seq} which is different than "
                           f"the translated seq {translated_original_seq} of the original sequence")
