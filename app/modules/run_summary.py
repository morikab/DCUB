import json
import os
import typing


class RunSummary(object):
    def __init__(self):
        self._run_summary = {}

    def add_to_run_summary(self, key: str, value: typing.Any) -> None:
        if key in self._run_summary:
            raise KeyError(F"Key {key} already exists in run summary")
        self._run_summary[key] = value

    def put_in_run_summary(self, key: str, value: typing.Any) -> None:
        self._run_summary[key] = value

    def append_to_run_summary(self, key: str, value: typing.Any) -> None:
        entry = self._run_summary.get(key, [])
        entry.append(value)
        self.put_in_run_summary(key=key, value=entry)

    def delete_from_run_summary(self, key:str) -> None:
        self._run_summary.pop(key)

    def save_run_summary(self, output_directory: str) -> None:
        output_path = os.path.join(output_directory, "run_summary.json")
        with open(output_path, "w") as output_file:
            json.dump(self._run_summary, output_file)

    def reset(self) -> None:
        self._run_summary = {}

    def get(self) -> typing.Dict[str, typing.Any]:
        return self._run_summary
