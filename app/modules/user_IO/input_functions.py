import os

import operator
import typing

from Bio import SeqFeature
from Bio import SeqIO
from Bio.Seq import Seq
import codonbias as cb
import pandas as pd

from logger_factory.logger_factory import LoggerFactory
from modules import models
from modules.configuration import Configuration
from modules.ORF.calculating_cai import relative_adaptiveness
from modules.ORF.TAI import TAI


logger = LoggerFactory.get_logger()
config = Configuration.get_config()


def tai_from_tgcnDB(org_name):
    base_path = os.path.join(os.path.dirname(__file__))
    tgcn_df = pd.read_csv(os.path.join(base_path, 'filtered_tgcn.csv'), index_col=0)
    all_org_tgcn_dict = tgcn_df.T.to_dict('list')

    if org_name in all_org_tgcn_dict.keys():
        tgcn_dict = dict(zip(tgcn_df.columns, all_org_tgcn_dict[org_name]))
        tai_weights = TAI(tgcn_dict).index
        logger.info(f'tGCN values were found, tAI profile was calculated')
    else:
        tai_weights = {}
        logger.info(f'tGCN values were not found, tAI profile was not calculated')
    return tai_weights


def extract_mrna_expression_levels(expression_csv_fid: str) -> typing.Tuple[typing.List[str], typing.List[float]]:
    expression_df = pd.read_csv(expression_csv_fid)

    gene_name_to_mrna_level = {}
    for idx, pair in enumerate(zip(expression_df.gene.to_list(), expression_df.mRNA_level.to_list())):
        measured_gene_name, expression_level = pair
        try:
            gene_name_to_mrna_level[measured_gene_name.lower()] = float(expression_level)
        except:
            continue
    mrna_levels = list(gene_name_to_mrna_level.values())
    mrna_names = list(gene_name_to_mrna_level.keys())

    return mrna_names, mrna_levels


def extract_protein_abundance_levels(expression_csv_fid: str) -> typing.Tuple[typing.List[str], typing.List[float]]:
    expression_df = pd.read_json(expression_csv_fid)

    gene_name_to_expression_level = {}
    for idx, pair in enumerate(zip(expression_df.name.to_list(), expression_df.abundance.to_list())):
        measured_gene_name, expression_level = pair
        try:
            gene_name_to_expression_level[measured_gene_name.lower()] = float(expression_level)
        except:
            continue
    protein_levels = list(gene_name_to_expression_level.values())
    protein_names = list(gene_name_to_expression_level.keys())

    return protein_names, protein_levels


def extract_expression_levels(expression_csv_fid: str,
                              expression_csv_type: str) -> typing.Tuple[typing.List[str], typing.List[float]]:
    expression_csv_type_mapping = {
        "mrna_levels": extract_mrna_expression_levels,
        "protein_abundance": extract_protein_abundance_levels,
    }
    logger.info(f"Extracting expression levels from: {expression_csv_type} file.")
    if expression_csv_type not in expression_csv_type_mapping:
        raise KeyError(F"Missing support for expression csv type: {expression_csv_type}")

    return expression_csv_type_mapping[expression_csv_type](expression_csv_fid)


def is_known_position_type(position: typing.Type[SeqFeature.Position]) -> bool:
    return isinstance(position, SeqFeature.ExactPosition) or isinstance(position, SeqFeature.BeforePosition) or \
           isinstance(position, SeqFeature.AfterPosition)


def is_known_location_type(feature) -> bool:
    for location_part in feature.location.parts:
        is_known_start_index = is_known_position_type(location_part.start)
        is_known_end_index = is_known_position_type(location_part.end)
        is_known_location_part = is_known_start_index and is_known_end_index

        if not is_known_location_part:
            return False

    return True


def does_have_only_exact_positions(feature) -> bool:
    return all(
        isinstance(location_part.start, SeqFeature.ExactPosition) and
        isinstance(location_part.end, SeqFeature.ExactPosition) for location_part in feature.location.parts
    )


def extract_gene_name(feature) -> typing.Optional[str]:
    name = feature.qualifiers.get("gene") or feature.qualifiers.get("locus_tag")
    if name is None:
        return name
    return name[0]
    # return "-".join(list(name))


def extract_cds(cds_features: typing.List, sequence: Seq) -> typing.Sequence[models.Cds]:
    cds_list = []
    for feature in cds_features:
        gene_name = extract_gene_name(feature)
        if gene_name is None:
            continue

        # TODO - should we remove this?
        if any(x.gene_name == gene_name for x in cds_list):
            continue

        cds = str(feature.extract(sequence))
        if len(cds) % 3 != 0:
            continue

        gene_function = " ".join(feature.qualifiers.get("product", []))
        cds_list.append(models.Cds(
            gene_name=gene_name,
            function=gene_function,
            sequence=cds,
        ))

    return cds_list

def extract_gene_data(genbank_path: str):
    gb_records = SeqIO.parse(genbank_path, format="genbank")
    all_cds = []
    for record in gb_records:
        if any(not is_known_location_type(x) for x in record.features):
            raise RuntimeError(f"Unknown location type found in {genbank_path}")

        genome = record.seq

        cds_features: typing.List[SeqFeature] = [
            x for x in record.features if x.type == "CDS" and does_have_only_exact_positions(x)
        ]
        all_cds.extend(extract_cds(cds_features, genome))
    return all_cds


def extract_gene_expression(
        cds: typing.Sequence[models.Cds],
        expression_csv_fid: typing.Optional[str] = None,
        expression_csv_type: typing.Optional[str] = None,
) -> typing.Optional[typing.Dict[str, float]]:
    should_use_expression_csv = expression_csv_fid is not None
    if not should_use_expression_csv:
        return None

    try:
        gene_expression_names, gene_expression_levels = extract_expression_levels(
            expression_csv_fid=expression_csv_fid,
            expression_csv_type=expression_csv_type,
        )
    except:
        logger.info("Expression data file is corrupt. \nMake sure that: ")
        logger.info("1. File is in csv format")
        logger.info("2. Gene names fit their NCBI naming ")
        logger.info("3. Column with the gene names is labeled 'gene' for mrna levels or 'name' for protein abundance ")
        logger.info("4. Column with the gene expression levels is labeled 'mRNA_level' or 'abundance', accordingly")
        return None

    estimated_expression = {}
    for cds_record in cds:
        try:
            index = gene_expression_names.index(cds_record.gene_name.lower())
            expression_level = gene_expression_levels[index]
            estimated_expression[cds_record.name_and_function] = expression_level
        except ValueError:
            continue
    return estimated_expression


def get_reference_genes_for_cai(
        cds_dict: typing.Dict[str, typing.Any],
        estimated_expression_dict: typing.Optional[typing.Dict[str, float]] = None,
) -> typing.Dict[str, str]:
    """
    calculates the cai weights - if estimated_expression dictionary has more than 3 times the number of ribosomal genes,
    15% most highly expressed genes will be used as reference set.
    in any other case, ribosomal genes will be used
    """
    ribosomal_proteins_count_threshold = config["INPUT"]["RIBOSOMAL_PROTEINS_COUNT_THRESHOLD"]
    ribosomal_proteins = {description: cds for description, cds in cds_dict.items() if "ribosom" in description}
    ribosomal_proteins_count = len(ribosomal_proteins)
    logger.info(F"Found {ribosomal_proteins_count} ribosomal proteins in input genome.")
    estimated_expression_dict = estimated_expression_dict or {}
    if len(estimated_expression_dict) < max(ribosomal_proteins_count, ribosomal_proteins_count_threshold) * 3:
        logger.info("Estimated expression dictionary does not have enough expression levels. CAI will be calculated "
                    "from a reference set of ribosomal proteins or the entire genome.")

        if ribosomal_proteins_count < ribosomal_proteins_count_threshold:
            logger.warning(F"Less than {ribosomal_proteins_count_threshold} ribosomal genes were found. This "
                           F"annotation is likely to be of low quality and therefore results may be less accurate."
                           F"CAI will be calculated from the entire genome as a reference set.")
            return cds_dict

        logger.info("CAI will be calculated from a reference set of ribosomal proteins.")
        return ribosomal_proteins

    logger.info("CAI will be calculated from a reference set of estimated expression dictionary.")
    logger.info(F"Expression levels were found for {len(estimated_expression_dict)}")
    estimated_expression_threshold = config["INPUT"]["EXPRESSION_PERCENTAGE_THRESHOLD"]
    sorted_estimated_expression = dict(
        sorted(estimated_expression_dict.items(), key=operator.itemgetter(1), reverse=True)
    )
    highly_expressed_genes_count = round(len(sorted_estimated_expression) * estimated_expression_threshold)
    logger.info(F"Calculate CAI weights from a reference set of {highly_expressed_genes_count} highly expressed "
                F"genes from estimated expression dictionary.")

    highly_expressed_names = list(sorted_estimated_expression.keys())[:highly_expressed_genes_count]
    reference_genes = {description: cds for description, cds in cds_dict.items() if description in
                       highly_expressed_names}
    return reference_genes


def calculate_tai_weights(organism_name: str) -> typing.Optional[cb.scores.TrnaAdaptationIndex]:
    # TODO - move to json file + pre load for multiple organisms or use the API to derive
    #  taxonomy level from the .gb file
    organism_name_to_url_mapping = {
        "Escherichia coli": "http://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/",
        "Bacillus subtilis": "http://gtrnadb.ucsc.edu/genomes/bacteria/Baci_subt_subtilis_168/"
    }
    if organism_name not in organism_name_to_url_mapping:
        logger.info(f"tGCN values were not found for {organism_name}, tAI profile was not calculated.")
        return None
    # TODO - consider using https://github.com/AliYoussef96/gtAI for better results
    logger.info(f"tGCN values were found for {organism_name}. Calculating tAI profile.")

    # Patch for fixing a certificate expiration issue when accessing gtrnadb
    import ssl
    ssl._create_default_https_context = ssl._create_unverified_context

    return cb.scores.TrnaAdaptationIndex(url=organism_name_to_url_mapping[organism_name], prokaryote=True)
