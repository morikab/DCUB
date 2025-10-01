import os
import random
import string
import typing
from pathlib import Path
from os import listdir
from os.path import isfile, join


current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.parent.resolve(), "example_data")

DEFAULT_MICROBIOME_PATH = os.path.join(base_path, 'arabidopsis_microbiome')
DEFAULT_ORGANISM_PRIORITY = 50
DEFAULT_CLUSTERS_COUNT = 1
DEFAULT_TUNING_PARAM = 0.5
DEFAULT_SEQUENCE_FILE_PATH = os.path.join(base_path, "mCherry_original.fasta")

DEFAULT_ECOLI_MRNA_EXPRESSION_LEVELS_FILE = "ecoli_mrna_level.csv"
DEFAULT_BACILLUS_MRNA_EXPRESSION_LEVELS_FILE = "bacillus_mrna_level.csv"

DEFAULT_ECOLI_PROTEIN_ABUNDANCE_FILE = "PA_ecoli.json"
DEFAULT_BACILLUS_PROTEIN_ABUNDANCE_FILE = "PA_bacillus.json"


def generate_random_string(length: int) -> str:
    letters_and_digits = string.ascii_letters + string.digits
    return "".join(random.choice(letters_and_digits) for _ in range(length))


def get_organisms_for_testing(
        organisms_count: int,
        percent_optimized: float = 0.5,
        microbiome_genome_path: str = DEFAULT_MICROBIOME_PATH,
        excluded_genomes: typing.Sequence[str] = None,
) -> typing.Tuple[typing.Sequence[str], typing.Sequence[str]]:
    wanted_hosts_count = round(organisms_count * percent_optimized)

    excluded_genomes = set(excluded_genomes) or ()
    genome_list = set(f for f in listdir(microbiome_genome_path) if isfile(join(microbiome_genome_path, f)))
    genome_list = genome_list.difference(excluded_genomes)

    if organisms_count > len(genome_list):
        raise ValueError("Not enough genomes in data. Select less genomes")

    selected_genomes = random.sample(genome_list, organisms_count)

    wanted_hosts = selected_genomes[:wanted_hosts_count]
    unwanted_hosts = selected_genomes[wanted_hosts_count:]

    return wanted_hosts, unwanted_hosts


def generate_testing_data(
        optimization_method: str,
        optimization_cub_index: str,
        wanted_hosts: typing.Sequence[str],
        unwanted_hosts: typing.Sequence[str],
        genome_path: str = DEFAULT_MICROBIOME_PATH,
        clusters_count: int = DEFAULT_CLUSTERS_COUNT,
        wanted_hosts_weights: typing.Dict[str, float] = None,
        unwanted_hosts_weights: typing.Dict[str, float] = None,
        tuning_param: float = DEFAULT_TUNING_PARAM,
        sequence_file_path: str = None,
        sequence: str = None,
        output_path: str = None,
) -> typing.Dict[str, typing.Any]:
    assert (sequence is not None or sequence_file_path is not None), \
        "Should provide either a sequence or a sequence file path"

    output_path = os.path.join("results", output_path)
    output_directory = os.path.join(output_path, F"{optimization_cub_index}_{optimization_method}_"
                                                 F"{len(wanted_hosts) + len(unwanted_hosts)}_"
                                                 F"{generate_random_string(4)}")
    # Path(output_directory).mkdir(parents=True, exist_ok=True)
    input_dict = {
        "sequence_file_path": sequence_file_path,
        "sequence": sequence,
        "tuning_param": tuning_param,
        "clusters_count": clusters_count,
        "optimization_method": optimization_method,
        "optimization_cub_index": optimization_cub_index,
        "output_path": output_directory,
        "organisms": {},
    }

    wanted_hosts_weights = wanted_hosts_weights or {}
    unwanted_hosts_weights = unwanted_hosts_weights or {}

    for host in wanted_hosts:
        input_dict["organisms"][host[:-3]] = {
            "genome_path": os.path.join(genome_path, host),
            "optimized": True,
            "expression_csv": None,
            "optimization_priority": wanted_hosts_weights.get(host) or DEFAULT_ORGANISM_PRIORITY,
        }

    for host in unwanted_hosts:
        input_dict["organisms"][host[:-3]] = {
            "genome_path": os.path.join(genome_path, host),
            "optimized": False,
            "expression_csv": None,
            "optimization_priority": unwanted_hosts_weights.get(host) or DEFAULT_ORGANISM_PRIORITY,
        }

    return input_dict


def generate_testing_data_for_ecoli_and_bacillus(
        optimization_method: str,
        optimization_cub_index: str,
        clusters_count: int = DEFAULT_CLUSTERS_COUNT,
        tuning_param: float = DEFAULT_TUNING_PARAM,
        is_ecoli_optimized: bool = False,
        sequence_file_path: str = None,
        sequence: str = None,
        output_path: str = None,
        should_use_mrna_levels: bool = False,
        should_use_protein_abundance: bool = True,
        evaluation_score: str = None,
):
    assert (sequence is not None or sequence_file_path is not None), \
        "Should provide either a sequence or a sequence file path"

    output_path = output_path or "results"
    output_directory = os.path.join(output_path, F"{optimization_cub_index}_{optimization_method}_ecoli_opt_"
                                                 F"{is_ecoli_optimized}_{generate_random_string(4)}")
    # Path(output_directory).mkdir(parents=True, exist_ok=True)

    if is_ecoli_optimized:
        opt_genome = "Escherichia-coli.gb"
        opt_mrna_levels = DEFAULT_ECOLI_MRNA_EXPRESSION_LEVELS_FILE if should_use_mrna_levels else None
        opt_protein_abundance = DEFAULT_ECOLI_PROTEIN_ABUNDANCE_FILE if should_use_protein_abundance else None
        deopt_genome = "Bacillus-subtilis.gb"
        deopt_mrna_levels = DEFAULT_BACILLUS_MRNA_EXPRESSION_LEVELS_FILE if should_use_mrna_levels else None
        deopt_protein_abundance = DEFAULT_BACILLUS_PROTEIN_ABUNDANCE_FILE if should_use_protein_abundance else None
    else:
        opt_genome = "Bacillus-subtilis.gb"
        opt_mrna_levels = DEFAULT_BACILLUS_MRNA_EXPRESSION_LEVELS_FILE if should_use_mrna_levels else None
        opt_protein_abundance = DEFAULT_BACILLUS_PROTEIN_ABUNDANCE_FILE if should_use_protein_abundance else None
        deopt_genome = "Escherichia-coli.gb"
        deopt_mrna_levels = DEFAULT_ECOLI_MRNA_EXPRESSION_LEVELS_FILE if should_use_mrna_levels else None
        deopt_protein_abundance = DEFAULT_ECOLI_PROTEIN_ABUNDANCE_FILE if should_use_protein_abundance else None

    inp_dict = {
        "sequence_file_path": sequence_file_path,
        "sequence": sequence,
        "tuning_param": tuning_param,
        "organisms": {},
        "clusters_count": clusters_count,
        "optimization_method": optimization_method,
        "optimization_cub_index": optimization_cub_index,
        "output_path": output_directory,
        "evaluation_score": evaluation_score or "average_distance",
    }

    expression_csv_type = None
    opt_expression_csv = None
    deopt_expression_csv = None
    if should_use_protein_abundance:
        expression_csv_type = "protein_abundance"
        opt_expression_csv = os.path.join(base_path, opt_protein_abundance)
        deopt_expression_csv = os.path.join(base_path, deopt_protein_abundance)
    elif should_use_mrna_levels:
        expression_csv_type = "mrna_levels"
        opt_expression_csv = os.path.join(base_path, opt_mrna_levels)
        deopt_expression_csv = os.path.join(base_path, deopt_mrna_levels)

    inp_dict['organisms'][opt_genome[:-3]] = {
        "genome_path": os.path.join(base_path, opt_genome),
        "optimized": True,
        "expression_csv_type": expression_csv_type,
        "expression_csv": opt_expression_csv,
        "optimization_priority": DEFAULT_ORGANISM_PRIORITY,
    }

    inp_dict["organisms"][deopt_genome[:-3]] = {
        "genome_path": os.path.join(base_path, deopt_genome),
        "optimized": False,
        "expression_csv_type": expression_csv_type,
        "expression_csv": deopt_expression_csv,
        "optimization_priority": DEFAULT_ORGANISM_PRIORITY,
    }
    return inp_dict
