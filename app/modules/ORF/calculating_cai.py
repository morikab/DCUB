from itertools import chain
from scipy.stats import gmean
from collections import Counter
# get rid of Biopython warning
import warnings
from Bio import BiopythonWarning

# TODO - consider using CAI code: from Bio.SeqUtils import CodonAdaptationIndex


warnings.simplefilter("ignore", BiopythonWarning)

genetic_code_dict = {
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


def _synonymous_codons(genetic_code_dict=genetic_code_dict):

    # invert the genetic code dictionary to map each amino acid to its codons
    codons_for_amino_acid = {}
    for codon, amino_acid in genetic_code_dict.items():
        codons_for_amino_acid[amino_acid] = codons_for_amino_acid.get(amino_acid, [])
        codons_for_amino_acid[amino_acid].append(codon)

    # create dictionary of synonymous codons
    # Example: {'CTT': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'ATG': ['ATG']...}
    return {
        codon: codons_for_amino_acid[genetic_code_dict[codon]]
        for codon in genetic_code_dict.keys()
    }


v = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y',
     'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
     'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
     'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
     'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
     'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D',
     'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
     'TGA': '-', 'TAA': '-', 'TAG': '-'}

_synonymous_codons = {
    'genetic_code': _synonymous_codons(v)}

_non_synonymous_codons = {
    k: {codon for codon in v.keys() if len(v[codon]) == 1}
    for k, v in _synonymous_codons.items()
}


def RSCU(sequences, genetic_code='genetic_code'):
    r"""Calculates the relative synonymous codon usage (RSCU) for a set of sequences.
    RSCU is 'the observed frequency of [a] codon divided by the frequency
    expected under the assumption of equal usage of the synonymous codons for an
    amino acid' (page 1283).
    In math terms, it is
    .. math::
        \frac{X_{ij}}{\frac{1}{n_i}\sum_{j=1}^{n_i}x_{ij}}
    "where :math:`X` is the number of occurrences of the :math:`j` th codon for
    the :math:`i` th amino acid, and :math:`n` is the number (from one to six)
    of alternative codons for the :math:`i` th amino acid" (page 1283).
    Args:
        sequences (list): The reference set of sequences.
        genetic_code (int, optional): The translation table to use. Defaults to 11, the standard genetic code.
    Returns:
        dict: The relative synonymous codon usage.
    Raises:
        ValueError: When an invalid sequence is provided or a list is not provided.
    """
    # ensure all input sequences are divisible by three
    sequences_new = []
    for sequence in sequences:
        if len(sequence) % 3 == 0:
            sequences_new.append(sequence)
        else:
            print(sequence)
        if not sequence:
            raise ValueError("Input sequence cannot be empty")
    sequences = sequences_new

    # count the number of each codon in the sequences
    sequences = (
        (sequence[i: i + 3].upper() for i in range(0, len(sequence), 3))
        for sequence in sequences
    )
    codons = chain.from_iterable(
        sequences
    )  # flat list of all codons (to be used for counting)
    counts = Counter(codons)

    # "if a certain codon is never used in the reference set... assign [its
    # count] a value of 0.5" (page 1285)
    for codon in v:
        if counts[codon] == 0:
            counts[codon] = 0.5

    # determine the synonymous codons for the genetic code
    synonymous_codons = _synonymous_codons[genetic_code]

    # hold the result as it is being calulated
    result = {}

    # calculate RSCU values
    for codon in v:
        result[codon] = counts[codon] / (
            (len(synonymous_codons[codon]) ** -1)
            * (sum((counts[_codon] for _codon in synonymous_codons[codon])))
        )

    return result


def relative_adaptiveness(sequences=None, RSCUs=None, genetic_code='genetic_code'):
    r"""Calculates the relative adaptiveness/weight of codons.
    The relative adaptiveness is "the frequency of use of that codon compared to
    the frequency of the optimal codon for that amino acid" (page 1283).
    In math terms, :math:`w_{ij}`, the weight for the :math:`j` th codon for
    the :math:`i` th amino acid is
    .. math::
        w_{ij} = \frac{\text{RSCU}_{ij}}{\text{RSCU}_{imax}}
    where ":math:`\text{RSCU}_{imax}` [is] the RSCU... for the frequently used
    codon for the :math:`i` th amino acid" (page 1283).
    Args:
        sequences (list, optional): The reference set of sequences.
        RSCUs (dict, optional): The RSCU of the reference set.
        genentic_code (int, optional): The translation table to use. Defaults to 11, the standard genetic code.
    Note:
        Either ``sequences`` or ``RSCUs`` is required.
    Returns:
        dict: A mapping between each codon and its weight/relative adaptiveness.
    Raises:
        ValueError: When neither ``sequences`` nor ``RSCUs`` is provided.
        ValueError: See :func:`RSCU` for details.
    """

    # ensure user gave only and only one input
    if sum([bool(sequences), bool(RSCUs)]) != 1:
        raise TypeError("Must provide either reference sequences or RSCU dictionary")

    # calculate the RSCUs if only given sequences
    if sequences:
        RSCUs = RSCU(sequences, genetic_code=genetic_code)

    # determine the synonymous codons for the genetic code
    synonymous_codons = _synonymous_codons[genetic_code]

    # calculate the weights
    weights = {}
    for codon in RSCUs:
        weights[codon] = RSCUs[codon] / max(
            (RSCUs[_codon] for _codon in synonymous_codons[codon])
        )

    return weights


def general_geomean(sequence_lst, weights, genetic_code='genetic_code'):
    scores = []
    for sequence in sequence_lst:
        sequence = sequence.upper()
        sequence = [sequence[i: i + 3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]
        sequence_weights = []
        for codon in sequence:
            if codon not in _non_synonymous_codons[genetic_code]:
                try:
                    if weights[codon] != 0:
                        sequence_weights.append(weights[codon])
                    else:
                        ValueError()
                except:
                    # if codon not in table (like if it conatians N or other ambigous chars) - it will be ignored
                    # (not counted in the seq length)
                    sequence_weights.append(sum(weights.values()) / len(weights.values()))
        scores.append(float(gmean(sequence_weights)))
    return scores
