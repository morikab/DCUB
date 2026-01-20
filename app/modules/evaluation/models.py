import typing
from dataclasses import dataclass

from modules import models


@dataclass
class EvaluationModuleResult:
    # TODO - Add additional stats (% optimized, % deoptimized, etc.)
    sequence: str
    average_distance_score: float
    weakest_link_score: float
    ratio_score: float

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "final_sequence": self.sequence,
            "average_distance_score": self.average_distance_score,
            "weakest_link_score": self.weakest_link_score,
            "ratio_score": self.ratio_score,
        }

    def get_score(self, score_type: models.EvaluationScore) -> float:
        score_value = score_type.value
        if score_value == models.EvaluationScore.average_distance.value:
            return self.average_distance_score
        if score_value == models.EvaluationScore.weakest_link.value:
            return self.weakest_link_score
        if score_value == models.EvaluationScore.ratio.value:
            return self.ratio_score
        raise ValueError(F"score type {score_type} (value: {score_value}) is not supported")
