from collections import defaultdict

from modules.run_summary import RunSummary
from modules.user_IO.input_functions import *
from modules.shared_functions_and_vars import DEFAULT_ORGANISM_PRIORITY
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons
from modules.shared_functions_and_vars import write_fasta


logger = LoggerFactory.get_logger()


class UserInputModule(object):
    def __init__(self, user_input: models.UserInput, skipped_codons_num: int):
        self.user_input = user_input
        self.skipped_codons_num = skipped_codons_num

    def run_module(
            self,
            run_summary: RunSummary,
    ) -> models.ModuleInput:
        logger.info('\n##########################')
        logger.info('# User Input #')
        logger.info('##########################')
        return self._parse_input(
            run_summary=run_summary,
        )

    def _parse_input(
            self,
            run_summary: RunSummary,
    ) -> models.ModuleInput:
        orf_sequence = self._parse_orf_sequence()

        organisms_list = self._parse_organisms_list(
            organisms_input_list=self.user_input.organisms,
            optimization_cub_index=self.user_input.orf_optimization_cub_index,
        )

        module_input = models.ModuleInput(
            organisms=organisms_list,
            sequence=orf_sequence,
            tuning_parameter=self.user_input.tuning_param,
            orf_optimization_method=self.user_input.orf_optimization_method,
            orf_optimization_cub_index=self.user_input.orf_optimization_cub_index,
            initiation_optimization_method=self.user_input.initiation_optimization_method,
            evaluation_score=self.user_input.evaluation_score_type,
            clusters_count=self.user_input.clusters_count,
            output_path=self.user_input.output_path,
            enable_hotspot_avoidance=self.user_input.enable_hotspot_avoidance,
            enable_motif_detection=self.user_input.enable_motif_detection,
        )

        run_summary.add_to_run_summary("module_input", module_input.summary)

        return module_input

    def _parse_orf_sequence(self) -> str:
        sequence = self.user_input.sequence
        if sequence is not None:
            return sequence

        orf_fasta_file_path = self.user_input.sequence_file_path
        logger.info(F"Sequence to be optimized given in the following file: {orf_fasta_file_path}")

        try:
            return str(SeqIO.read(orf_fasta_file_path, "fasta").seq)
        except:
            raise ValueError(
                F"Error in orf sequence .fasta file: {orf_fasta_file_path}. Make sure you inserted a valid .fasta file "
                F"containing a single record")

    def _parse_organisms_list(
            self,
            organisms_input_list: typing.Dict[str, models.OrganismRequest],
            optimization_cub_index: models.ORFOptimizationCubIndex,
    ) -> typing.List[models.Organism]:
        organisms_list = []
        organisms_names = set()
        for organism_key, organism_input in organisms_input_list.items():
            try:
                organism = self._parse_single_organism_input(
                    organism_input=organism_input,
                    optimization_cub_index=optimization_cub_index,
                )
            except Exception as e:
                raise ValueError(f"Error in organism input: {organism_key}, re-check your input")
            if organism.name in organisms_names:
                raise ValueError(f"Organism: {organism.name}'s genome is inserted twice, re-check your input.")
            organisms_list.append(organism)
            organisms_names.add(organism.name)

        # Normalize prioritization weights - from user's defined values (1-100) to normalized values in range (0, 1)
        logger.info("-----------------------------")
        logger.info("Normalized Prioritization Weights")
        logger.info("-----------------------------")
        total_optimized_weights = sum([organism.optimization_priority for organism in organisms_list if
                                       organism.is_optimized])
        total_deoptimized_weights = sum([organism.optimization_priority for organism in organisms_list if
                                         not organism.is_optimized])
        for organism in organisms_list:
            if organism.is_optimized:
                organism.optimization_priority /= total_optimized_weights
            else:
                organism.optimization_priority /= total_deoptimized_weights
            logger.info(f"{organism.name} : {organism.optimization_priority}")

        return organisms_list

    def _parse_single_organism_input(
            self,
            organism_input: models.OrganismRequest,
            optimization_cub_index: models.ORFOptimizationCubIndex,
    ) -> models.Organism:

        is_optimized = organism_input.optimized
        gb_path = organism_input.genome_path

        # FIXME - delete
        parsed_organism_file_name = f"{gb_path.strip('.gb')}_{is_optimized}_parsed"
        parsed_organism_file = parsed_organism_file_name + ".json"

        # FIXME - end
        try:
            with open(gb_path, "r") as gb_file_handle:
                first_gb_record = next(SeqIO.parse(gb_file_handle, format="genbank"))
            organism_name = " ".join(first_gb_record.description.split()[:2])
        except Exception as e:
            raise ValueError(
                f'Error in genome GenBank file: {gb_path}, make sure you inserted an undamaged .gb file containing '
                f'the full genome sequence and annotations'
            ) from e
        logger.info("------------------------------------------")
        logger.info(f"Parsing information for {organism_name}:")
        logger.info(f"Raw organism request is: {organism_input}")
        logger.info("------------------------------------------")
        logger.info(f"Organism is defined as {'wanted' if is_optimized else 'unwanted'}")
        cds = extract_gene_data(genbank_path=gb_path)

        estimated_expression = extract_gene_expression(
            cds=cds,
            expression_file_path=organism_input.expression_file_path,
            expression_data_type=organism_input.expression_data_type,
            expression_file_format=organism_input.expression_file_format,
        )

        cds_dict = {cds_record.name_and_function: cds_record.sequence for cds_record in cds}
        gene_names = list(cds_dict.keys())

        cai_weights = None
        cai_scores_dict = None
        tai_weights = None
        tai_scores_dict = None
        reference_genes = None

        if optimization_cub_index.is_codon_adaptation_index:
            reference_genes = get_reference_genes_for_cai(cds_dict, estimated_expression)
            reference_genes_for_cub_calculation = [seq[self.skipped_codons_num*3:] for seq in reference_genes.values()]
            cai = cb.scores.CodonAdaptationIndex(ref_seq=reference_genes_for_cub_calculation, ignore_stop=False)
            cai_weights = cai.weights.to_dict()
            cai_scores_dict = {gene_name: cai.get_score(cds_dict[gene_name]) for gene_name in gene_names}

        if optimization_cub_index.is_trna_adaptation_index:
            tai = calculate_tai_weights(organism_name)
            if tai is not None:
                tai_weights = tai.weights.to_dict()
                tai_scores_dict = {gene_name: tai.get_score(cds_dict[gene_name][self.skipped_codons_num*3:]) for gene_name in gene_names if len(cds_dict[gene_name]) > (self.skipped_codons_num+1)*3}

        codon_frequencies = _calculate_codon_frequencies(cds_dict)

        optimization_priority = organism_input.optimization_priority or DEFAULT_ORGANISM_PRIORITY
        organism_object = models.Organism(name=organism_name,
                                          cai_profile=cai_weights,
                                          tai_profile=tai_weights,
                                          codon_frequencies=codon_frequencies,
                                          cai_scores=cai_scores_dict,
                                          tai_scores=tai_scores_dict,
                                          reference_genes=list(reference_genes.keys()) if reference_genes else None,
                                          is_optimized=is_optimized,
                                          optimization_priority=optimization_priority)

        # FIXME - delete
        org_summary = organism_object.summary
        org_summary["cai_scores"] = cai_scores_dict
        org_summary["tai_scores"] = tai_scores_dict
        org_summary["cds_dict"] = cds_dict
        org_summary["reference_genes"] = reference_genes
        write_fasta(fid=organism_name, list_seq=list(cds_dict.values()), list_name=list(cds_dict.keys()))

        # with open(parsed_organism_file, "w") as organism_file:
        #     json.dump(org_summary, organism_file)
        # FIXME - end

        if optimization_cub_index.is_codon_adaptation_index:
            logger.info(F"cai_std={organism_object.cai_std}, cai_avg={organism_object.cai_avg}")
        if optimization_cub_index.is_trna_adaptation_index:
            logger.info(F"tai_std={organism_object.tai_std}, tai_avg={organism_object.tai_avg}")
        return organism_object


def _calculate_codon_frequencies(reference_cds):
    codons_counter = defaultdict(int)

    for cds in reference_cds.values():
        if len(cds) % 3 != 0:
            continue
        for i in range(0, len(cds), 3):
            codon = cds[i:i + 3]
            codons_counter[codon] += 1

    max_aa_freq = {}
    # Calculate relative frequncies
    for amino_acid, codons in synonymous_codons.items():
        total_amino_acid_codons = 0
        for amino_acid_codon in codons:
            total_amino_acid_codons += codons_counter[amino_acid_codon]

        # In case a certain amino acid is missing from the reference cds collection, no need to normalize
        if total_amino_acid_codons == 0:
            continue

        for amino_acid_codon in codons:
            codons_counter[amino_acid_codon] /= total_amino_acid_codons

    return codons_counter
