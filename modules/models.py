import statistics
import typing
from dataclasses import dataclass
from enum import Enum


@dataclass
class Cds:
    gene_name: str
    function: str
    sequence: str

    @property
    def name_and_function(self) -> str:
        return f"{self.gene_name}|{self.function}"


class Organism(object):
    CAI_PROFILE_ATTRIBUTE_NAME = "cai_profile"
    TAI_PROFILE_ATTRIBUTE_NAME = "tai_profile"

    def __init__(self,
                 name: str,
                 is_optimized: bool,
                 optimization_priority: float,
                 cai_profile: typing.Optional[typing.Dict] = None,
                 tai_profile: typing.Optional[typing.Dict] = None,
                 codon_frequencies: typing.Optional[typing.Dict] = None,
                 cai_scores: typing.Optional[typing.Dict] = None,
                 tai_scores: typing.Optional[typing.Dict] = None,
                 reference_genes: typing.Optional[typing.Sequence] = None):
        self.name = name
        self.cai_profile = cai_profile
        self.tai_profile = tai_profile
        self.codon_frequencies = codon_frequencies
        self.cai_scores = cai_scores
        self.tai_scores = tai_scores
        self.reference_genes = reference_genes
        self.is_optimized = is_optimized
        self.optimization_priority = optimization_priority

    def filter_reference_genes(self, scores: typing.Dict[str, float]):
        return [scores[gene_name] for gene_name in scores.keys() if gene_name not in self.reference_genes]

    @property
    def cai_avg(self) -> typing.Optional[float]:
        if self.cai_scores is None:
            return None
        scores = self.cai_scores
        return statistics.mean(scores.values())

    @property
    def tai_avg(self) -> float:
        return statistics.mean(self.tai_scores.values()) if self.tai_scores else None

    @property
    def cai_std(self) -> typing.Optional[float]:
        if self.cai_scores is None:
            return None
        scores = self.cai_scores
        return statistics.stdev(scores.values())

    @property
    def tai_std(self) -> float:
        return statistics.stdev(self.tai_scores.values()) if self.tai_scores else None

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "name": self.name,
            "is_wanted": self.is_optimized,
            "optimization_priority": self.optimization_priority,
            "cai_weights": self.cai_profile,
            "codon_frequencies": self.codon_frequencies,
            "cai_avg": self.cai_avg,
            "cai_std": self.cai_std,
            "tai_weights": self.tai_profile,
            "tai_avg": self.tai_avg,
            "tai_std": self.tai_std,
        }


class ORFOptimizationMethod(Enum):
    single_codon_ratio = "single_codon_ratio"
    single_codon_diff = "single_codon_diff"
    single_codon_weakest_link = "single_codon_weakest_link"
    zscore_single_aa_ratio = "zscore_single_aa_ratio"
    zscore_bulk_aa_ratio = "zscore_bulk_aa_ratio"
    zscore_single_aa_diff = "zscore_single_aa_diff"
    zscore_bulk_aa_diff = "zscore_bulk_aa_diff"
    zscore_single_aa_weakest_link = "zscore_single_aa_weakest_link"
    zscore_bulk_aa_weakest_link = "zscore_bulk_aa_weakest_link"
    single_organism = "single_wanted_organism"

    @property
    def is_single_codon_optimization(self) -> bool:
        return self in (ORFOptimizationMethod.single_codon_ratio,
                        ORFOptimizationMethod.single_codon_diff,
                        ORFOptimizationMethod.single_codon_weakest_link)

    @property
    def is_zscore_single_aa_optimization(self) -> bool:
        return self in (ORFOptimizationMethod.zscore_single_aa_ratio,
                        ORFOptimizationMethod.zscore_single_aa_diff,
                        ORFOptimizationMethod.zscore_single_aa_weakest_link)

    @property
    def is_zscore_bulk_aa_optimization(self) -> bool:
        return self in (ORFOptimizationMethod.zscore_bulk_aa_ratio,
                        ORFOptimizationMethod.zscore_bulk_aa_diff,
                        ORFOptimizationMethod.zscore_bulk_aa_weakest_link)

    @property
    def is_zscore_optimization(self) -> bool:
        return self.is_zscore_single_aa_optimization or self.is_zscore_bulk_aa_optimization

    @property
    def is_zscore_ratio_score_optimization(self) -> bool:
        return self in (ORFOptimizationMethod.zscore_single_aa_ratio, ORFOptimizationMethod.zscore_bulk_aa_ratio)

    @property
    def is_zscore_weakest_link_score_optimization(self) -> bool:
        return self in (ORFOptimizationMethod.zscore_single_aa_weakest_link,
                        ORFOptimizationMethod.zscore_bulk_aa_weakest_link)

    @property
    def is_zscore_diff_score_optimization(self) -> bool:
        return self in (ORFOptimizationMethod.zscore_single_aa_diff, ORFOptimizationMethod.zscore_bulk_aa_diff)


class InitiationOptimizationMethod(Enum):
    original = "original"
    external_module = "external"
    weak_folding = "weak_folding"


class ORFOptimizationCubIndex(Enum):
    codon_adaptation_index = "CAI"
    trna_adaptation_index = "tAI"
    max_codon_trna_adaptation_index = "max_CAI_tAI"

    @property
    def is_codon_adaptation_index(self) -> bool:
        return self in (ORFOptimizationCubIndex.codon_adaptation_index,
                        ORFOptimizationCubIndex.max_codon_trna_adaptation_index)

    @property
    def is_trna_adaptation_index(self) -> bool:
        return self in (ORFOptimizationCubIndex.trna_adaptation_index,
                        ORFOptimizationCubIndex.max_codon_trna_adaptation_index)


class EvaluationScore(Enum):
    average_distance = "average_distance"
    weakest_link = "weakest_link"
    ratio = "ratio"


@dataclass
class ModuleInput:
    organisms: typing.List[Organism]
    sequence: str
    output_path: str
    tuning_parameter: float
    clusters_count: int
    orf_optimization_method: ORFOptimizationMethod = ORFOptimizationMethod.zscore_bulk_aa_ratio
    orf_optimization_cub_index: ORFOptimizationCubIndex = ORFOptimizationCubIndex.max_codon_trna_adaptation_index
    initiation_optimization_method: InitiationOptimizationMethod = InitiationOptimizationMethod.original
    evaluation_score: EvaluationScore = EvaluationScore.average_distance

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "sequence": self.sequence,
            "tuning_parameter": self.tuning_parameter,
            "orf_optimization_method": self.orf_optimization_method.value,
            "orf_optimization_cub_index": self.orf_optimization_cub_index.value,
            "initiation_optimization_method": self.initiation_optimization_method.value,
            "evaluation_score": self.evaluation_score.value,
            "organisms": [organism.summary for organism in self.organisms],
        }


@dataclass
class SequenceZscores:
    initial_wanted_hosts_scores: typing.List[float] = None
    wanted_hosts_weights: typing.List[float] = None
    initial_unwanted_hosts_scores: typing.List[float] = None
    unwanted_hosts_weights: typing.List[float] = None

    min_score_for_normalization: typing.Optional[float] = None
    max_score_for_normalization: typing.Optional[float] = None
    normalized_wanted_hosts_scores: typing.Optional[typing.List[float]] = None
    normalized_unwanted_hosts_scores: typing.Optional[typing.List[float]] = None

    def normalize(self, min_zscore: float, max_zscore: float) -> None:
        epsilon = 1e-6

        self.normalized_wanted_hosts_scores = [
            (score - min_zscore + epsilon) / (max_zscore - min_zscore + epsilon) for
            score in self.wanted_hosts_scores
        ]
        self.normalized_unwanted_hosts_scores = [
            (score - min_zscore + epsilon) / (max_zscore - min_zscore + epsilon) for
            score in self.unwanted_hosts_scores
        ]
        self.min_score_for_normalization = min_zscore
        self.max_score_for_normalization = max_zscore

    @property
    def min_zscore(self) -> float:
        return min(min(self.wanted_hosts_scores), min(self.unwanted_hosts_scores))

    @property
    def max_zscore(self) -> float:
        return max(max(self.wanted_hosts_scores), max(self.unwanted_hosts_scores))

    @property
    def wanted_hosts_scores(self) -> typing.List[float]:
        return self.normalized_wanted_hosts_scores or self.initial_wanted_hosts_scores

    @property
    def unwanted_hosts_scores(self) -> typing.List[float]:
        return self.normalized_unwanted_hosts_scores or self.initial_unwanted_hosts_scores

    @property
    def all_scores(self) -> list[float]:
        return self.wanted_hosts_scores + self.unwanted_hosts_scores
