import os
import yaml
from pathlib import Path


class Configuration(object):
    _CONFIG_PATH = str(os.path.join(Path(__file__).parent.resolve(), "configuration.yaml"))
    _config = None

    @classmethod
    def read_config(cls):
        with open(cls._CONFIG_PATH, "r") as config_file:
            return yaml.safe_load(config_file)

    @classmethod
    def get_config(cls):
        if cls._config is None:
            cls._config = cls.read_config()
        return cls._config
