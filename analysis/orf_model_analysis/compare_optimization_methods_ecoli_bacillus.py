import json
import os
import re
import typing
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from analysis.input_testing_data.generate_input_testing_data_for_modules import \
    generate_testing_data_for_ecoli_and_bacillus
from modules.main import run_modules
from modules.main import run_input_processing
from modules.main import run_orf_module


current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
DEFAULT_SEQUENCE_FILE_PATH = os.path.join(base_path, "mCherry_original.fasta")


def run_from_fasta_file(fasta_file_path: str,
                        records_file_path: str,
                        start_record: str,
                        max_records_count: int) -> None:
    with open(fasta_file_path, "r") as fasta_handle:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    with open(records_file_path, "r") as records_file:
        records = json.load(records_file)

    is_record_found = False
    count = 0
    for record in records.values():
        if record == start_record:
            is_record_found = True
        if not is_record_found:
            continue
        if count >= max_records_count:
            break

        count += 1
        value = genome_dict[record]
        run_all_methods(orf_sequence=str(value.seq),
                        output_path=record.replace('|', '-'))


def run_for_endogenous_genes(fasta_file_path: str,
                             optimization_method: str,
                             optimization_cub_index: str,
                             is_ecoli_optimized: bool,
                             output_path: str) -> None:
    with open(fasta_file_path, "r") as fasta_handle:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    results_dict = {}
    for gene_name, gene_sequence in genome_dict.items():
        gene_sequence = str(gene_sequence.seq)
        if len(gene_sequence) % 3 != 0:
            print(F"Invalid length {len(gene_sequence)} for gene {gene_name}")
            continue
        results_dict[gene_name] = run_single_method_for_orf_sequence(
            optimization_method=optimization_method,
            optimization_cub_index=optimization_cub_index,
            is_ecoli_optimized=is_ecoli_optimized,
            output_path=output_path,
            orf_sequence=gene_sequence)

    with open(
            F"{optimization_cub_index}_{optimization_method}_{Path(fasta_file_path).name[:10]}_{is_ecoli_optimized}_fasta_results.json",
            "w") as results_file:
        json.dump(results_dict, results_file)


def run_all_methods(orf_sequence: typing.Optional[str] = None,
                    orf_sequence_file: typing.Optional[str] = None,
                    output_path: typing.Optional[str] = None):
    for optimization_method in [
        "single_codon_ratio",
        "single_codon_diff",
        "single_codon_weakest_link",
        # "zscore_single_aa_ratio",
        "zscore_bulk_aa_ratio",
        # "zscore_single_aa_diff",
        "zscore_bulk_aa_diff",
        # "zscore_single_aa_weakest_link",
        "zscore_bulk_aa_weakest_link",
    ]:
        for optimization_cub_index in ["CAI", "tAI"]:
            for direction in [True, False]:
                run_single_method_for_orf_sequence(optimization_method=optimization_method,
                                                   optimization_cub_index=optimization_cub_index,
                                                   is_ecoli_optimized=direction,
                                                   orf_sequence=orf_sequence,
                                                   orf_sequence_file=orf_sequence_file,
                                                   output_path=output_path)


def extract_ncbi_sequences_for_analysis(fasta_file_path: str) -> None:
    with open(fasta_file_path, "r") as fasta_handle:
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

    missing_genes = []
    gene_mapping = defaultdict(list)
    description_to_gene_mapping = {}

    for record, value in genome_dict.items():
        if "NC" not in record or "XP" in record:
            continue
        parameters = re.findall("gene=.*]", value.description)
        if not parameters:
            missing_genes.append(record)
            continue
        gene_parameter = parameters[0].split("]")[0]
        gene_name = gene_parameter.strip("gene=")

        if len(value.seq) % 3 != 0:
            print(F"Skip sequence {record} for {gene_name} because of length not divisible by 3.")
            continue

        gene_mapping[gene_name].append(record)
        description_to_gene_mapping[record] = gene_name

    with open("gene_mapping.json", "w") as genes_file:
        json.dump(gene_mapping, genes_file)
    with open("description_to_gene_mapping.json", "w") as description_to_genes_file:
        json.dump(description_to_gene_mapping, description_to_genes_file)
    with open("missing_genes.txt", "w") as missing_genes_file:
        for gene in missing_genes:
            missing_genes_file.write(gene)

    gene_to_longest_sequence = {
        key: max(value, key=lambda x: len(genome_dict[x].seq)) for key, value in gene_mapping.items()
    }
    with open("gene_to_longest_sequence.json", "w") as genes_file:
        json.dump(gene_to_longest_sequence, genes_file)


def run_single_method_for_orf_sequence(optimization_method: str,
                                       is_ecoli_optimized: bool,
                                       orf_sequence: typing.Optional[str] = None,
                                       orf_sequence_file: typing.Optional[str] = None,
                                       output_path: typing.Optional[str] = None,
                                       optimization_cub_index: str = "CAI",
                                       tuning_param: float = 0.5,
                                       initiation_optimization_method: str = "original"):
    default_user_inp_raw = generate_testing_data_for_ecoli_and_bacillus(
        orf_optimization_method=optimization_method,
        orf_optimization_cub_index=optimization_cub_index,
        clusters_count=1,
        tuning_param=tuning_param,
        is_ecoli_optimized=is_ecoli_optimized,
        sequence=orf_sequence,
        sequence_file_path=orf_sequence_file,
        output_path=os.path.join("results", output_path),
        initiation_optimization_method=initiation_optimization_method,  # original, external
    )
    return run_modules(default_user_inp_raw, should_run_output_module=False)
    # return run_orf_module(default_user_inp_raw)
    # run_input_processing(default_user_inp_raw)


def compare_gene_mappings() -> None:
    with open("gene_to_longest_sequence_all.json") as gene_mapping_file:
        previous_mapping = json.load(gene_mapping_file)

    with open("gene_to_longest_sequence.json") as gene_mapping_file:
        current_mapping = json.load(gene_mapping_file)

    changed_genes = []
    uncovered_genes = []
    for gene in previous_mapping.keys():
        if current_mapping.get(gene) is None:
            uncovered_genes.append(gene)
            continue

        if previous_mapping[gene] != current_mapping[gene]:
            changed_genes.append(current_mapping[gene])

    with open("changed_genes.txt", "w") as changed_genes_file:
        for gene in changed_genes:
            changed_genes_file.write(current_mapping[gene]+"\n")


def generate_sequences_fasta_file(root_dir) -> None:
    filename = "run_summary.json"

    sequences = []
    sequences_names = []
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file == filename:
                directory_name = Path(root).name
                file_path = os.path.join(root, file)
                with open(file_path, "r") as summary_file:
                    results_json = json.load(summary_file)

                seq = results_json["evaluation"]["final_sequence"]
                sequences.append(seq)
                sequences_names.append(directory_name[:-5])

    from modules.shared_functions_and_vars import write_fasta

    write_fasta(os.path.join(root_dir, "mcherry_variants"), sequences, sequences_names)


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Analysis script parser")
    # parser.add_argument('-s', '--start', type=str, required=True, help="Fasta record description to start running from")
    # parser.add_argument('-n', '--number', type=int, help="Number of records to parse from the given start record")
    # parser.add_argument('-f', '--fasta', type=str, help="Fasta file for orf sequences to run on (for endogenous run)")
    # parser.add_argument('-m', '--method', type=str, help="Optimization method")
    # parser.add_argument('-i', '--index', type=str, help="Optimization CUB index")
    # parser.add_argument('--opt', type=bool, help="Boolean indicating whether e.coli is optimized or not")
    # parser.add_argument('-o', '--output', type=str, help="Output path")
    #
    # args = parser.parse_args()

    # run_all_methods(orf_sequence_file=DEFAULT_SEQUENCE_FILE_PATH,
    #                 output_path="mcherry_debug")

    # run_single_method_for_orf_sequence(optimization_method="zscore_single_aa_ratio",
    #                                    optimization_cub_index="CAI",
    #                                    is_ecoli_optimized=True,
    #                                    orf_sequence_file=DEFAULT_SEQUENCE_FILE_PATH,
    #                                    output_path="mcherry_debug")

    # fasta_file_path = r"C:\projects\Igem_TAU_2021_moran\analysis\example_data\Escherichia-coli.fasta"
    # with open(fasta_file_path, "r") as fasta_handle:
    #     genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"), lambda r: r.description)

    # with open(r"C:\projects\Igem_TAU_2021_moran\analysis\example_data\Bacillus-subtilis.fasta", "r") as fasta_handle:
    #     genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"), lambda r: r.description)
    # gene_name = "thrA|fused aspartate kinase/homoserine dehydrogenase 1"
    #
    # gene_sequence = genome_dict[gene_name]
    # gene_sequence = str(gene_sequence.seq)
    initiation_type = "external"
    results = run_single_method_for_orf_sequence(
        # optimization_method="zscore_single_aa_diff",
        optimization_method="zscore_bulk_aa_diff",
        # optimization_method="single_codon_diff",
        optimization_cub_index="CAI",
        is_ecoli_optimized=False,
        output_path=f"mcherry-debug",
        orf_sequence_file=DEFAULT_SEQUENCE_FILE_PATH,
        initiation_optimization_method=initiation_type,
        # orf_sequence="",
    )
    # FIXME - debugging sequences 5 & 7 for in-vitro experiment
    print("Final seq: " + results["final_evaluation"]["final_sequence"])
    # TODO - add to the logs the optimal and second best codons and make the replacement manually
