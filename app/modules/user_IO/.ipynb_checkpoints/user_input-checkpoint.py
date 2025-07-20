import json
from collections import defaultdict

from modules.run_summary import RunSummary
from modules.user_IO.input_functions import *
from modules.shared_functions_and_vars import DEFAULT_ORGANISM_PRIORITY
from modules.shared_functions_and_vars import nt_to_aa
from modules.shared_functions_and_vars import synonymous_codons
from modules.shared_functions_and_vars import write_fasta


logger = LoggerFactory.get_logger()


class UserInputModule(object):
    @staticmethod
    def get_name() -> str:
        return "User Input"

    @classmethod
    def run_module(
            cls,
            user_inp_raw: typing.Dict,
            initiation_optimized_codons_num: int,
            run_summary: RunSummary,
    ) -> models.ModuleInput:
        logger.info('\n##########################')
        logger.info('# User Input #')
        logger.info('##########################')
        return cls._parse_input(
            module_input=user_inp_raw,
            skipped_codons_num=initiation_optimized_codons_num,
            run_summary=run_summary,
        )

    @classmethod
    def _parse_input(
            cls,
            module_input: typing.Dict[str, typing.Any],
            skipped_codons_num: int,
            run_summary: RunSummary,
    ) -> models.ModuleInput:
        orf_optimization_cub_index = models.ORFOptimizationCubIndex(module_input["orf_optimization_cub_index"]) if \
            module_input.get("orf_optimization_cub_index") else None
        orf_optimization_method = models.ORFOptimizationMethod(module_input["orf_optimization_method"]) if \
            module_input.get("orf_optimization_method") else None
        initiation_optimization_method = models.InitiationOptimizationMethod(
            module_input["initiation_optimization_method"]
        ) if module_input.get("initiation_optimization_method") else None
        evaluation_score = models.EvaluationScore(module_input["evaluation_score"]) if \
            module_input.get("evaluation_score") else None
        tuning_parameter = module_input["tuning_param"]
        clusters_count = module_input["clusters_count"]
        output_path = module_input.get("output_path")

        orf_sequence = cls._parse_orf_sequence(module_input)
        logger.info(F"Open reading frame sequence for optimization is:\n{orf_sequence}")

        organisms_list = cls._parse_organisms_list(organisms_input_list=module_input["organisms"],
                                                   optimization_cub_index=orf_optimization_cub_index,
                                                   skipped_codons_num=skipped_codons_num)

        module_input = models.ModuleInput(organisms=organisms_list,
                                          sequence=orf_sequence,
                                          tuning_parameter=tuning_parameter,
                                          orf_optimization_method=orf_optimization_method,
                                          orf_optimization_cub_index=orf_optimization_cub_index,
                                          initiation_optimization_method=initiation_optimization_method,
                                          evaluation_score=evaluation_score,
                                          clusters_count=clusters_count,
                                          output_path=output_path)

        run_summary.add_to_run_summary("module_input", module_input.summary)

        return module_input

    @classmethod
    def _parse_orf_sequence(cls, user_input: typing.Dict[str, typing.Any]) -> str:
        sequence = user_input.get("sequence")
        if sequence is not None:
            return sequence

        orf_fasta_file_path = user_input["sequence_file_path"]
        logger.info(F"Sequence to be optimized given in the following file: {orf_fasta_file_path}")

        try:
            return str(SeqIO.read(orf_fasta_file_path, "fasta").seq)
        except:
            raise ValueError(
                F"Error in orf sequence .fasta file: {orf_fasta_file_path}. Make sure you inserted a valid .fasta file "
                F"containing a single record")

    @classmethod
    def _parse_organisms_list(cls,
                              organisms_input_list: typing.Dict[str, typing.Any],
                              optimization_cub_index: models.ORFOptimizationCubIndex,
                              skipped_codons_num: int) -> typing.List[models.Organism]:
        organisms_list = []
        organisms_names = set()
        for organism_key, organism_input in organisms_input_list.items():
            try:
                organism = cls._parse_single_organism_input(
                    organism_input=organism_input,
                    optimization_cub_index=optimization_cub_index,
                    skipped_codons_num=skipped_codons_num,
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

    @staticmethod
    def _parse_single_organism_input(organism_input: typing.Dict[str, typing.Any],
                                     optimization_cub_index: models.ORFOptimizationCubIndex,
                                     skipped_codons_num: int) -> models.Organism:

        is_optimized = organism_input["optimized"]
        gb_path = organism_input["genome_path"]

        # FIXME - delete
        parsed_organism_file_name = f"{gb_path.strip('.gb')}_{is_optimized}_parsed"
        parsed_organism_file = parsed_organism_file_name + ".json"

        # FIXME - delete
        # if os.path.exists(parsed_organism_file):
        #     with open(parsed_organism_file) as org_file:
        #         organism_data = json.load(org_file)
        #         return models.Organism(name=organism_data["name"],
        #                                cai_profile=organism_data["cai_weights"],
        #                                tai_profile=organism_data["tai_weights"],
        #                                cai_scores=organism_data["cai_scores"],
        #                                tai_scores=organism_data["tai_scores"],
        #                                reference_genes=organism_data["reference_genes"],
        #                                is_optimized=organism_data["is_wanted"],
        #                                optimization_priority=organism_data["optimization_priority"])

        # FIXME - end
        try:
            gb_file = SeqIO.read(gb_path, format="gb")
            organism_name = " ".join(gb_file.description.split()[:2])
        except:
            raise ValueError(
                f'Error in genome GenBank file: {gb_path}, make sure you inserted an undamaged .gb file containing '
                f'the full genome sequence and annotations'
            )
        logger.info("------------------------------------------")
        logger.info(f"Parsing information for {organism_name}:")
        logger.info("------------------------------------------")
        logger.info(f"Organism is defined as {'wanted' if is_optimized else 'unwanted'}")
        cds = extract_gene_data(genbank_path=gb_path)

        exp_csv_type = organism_input['expression_csv_type']
        exp_csv_fid = organism_input['expression_csv']
        estimated_expression = extract_gene_expression(cds=cds,
                                                       expression_csv_fid=exp_csv_fid,
                                                       expression_csv_type=exp_csv_type)

        cds_dict = {cds_record.name_and_function: cds_record.sequence for cds_record in cds}
        gene_names = list(cds_dict.keys())

        cai_weights = None
        cai_scores_dict = None
        tai_weights = None
        tai_scores_dict = None
        reference_genes = None

        if optimization_cub_index.is_codon_adaptation_index:
            reference_genes = get_reference_genes_for_cai(cds_dict, estimated_expression)
            reference_genes_for_cub_calculation = [seq[skipped_codons_num*3:] for seq in reference_genes.values()]
            cai = cb.scores.CodonAdaptationIndex(ref_seq=reference_genes_for_cub_calculation, ignore_stop=False)
            cai_weights = cai.weights.to_dict()
            cai_scores_dict = {gene_name: cai.get_score(cds_dict[gene_name]) for gene_name in gene_names}

        if optimization_cub_index.is_trna_adaptation_index:
            tai = calculate_tai_weights(organism_name)
            if tai is not None:
                tai_weights = tai.weights.to_dict()
                tai_scores_dict = {gene_name: tai.get_score(cds_dict[gene_name][skipped_codons_num*3:]) for gene_name in gene_names if len(cds_dict[gene_name]) > (skipped_codons_num+1)*3}

        codon_frequencies = _calculate_codon_frequencies(cds_dict)

        optimization_priority = organism_input.get("optimization_priority") or DEFAULT_ORGANISM_PRIORITY
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

        with open(parsed_organism_file, "w") as organism_file:
            json.dump(org_summary, organism_file)
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
