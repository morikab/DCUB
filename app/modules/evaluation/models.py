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
        if score_type == models.EvaluationScore.average_distance:
            return self.average_distance_score
        if score_type == models.EvaluationScore.weakest_link:
            return self.weakest_link_score
        if score_type == models.EvaluationScore.ratio:
            return self.ratio_score
        raise ValueError(F"score type {score_type} is not supported")
