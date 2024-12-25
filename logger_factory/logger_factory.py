import logging
import os
import typing
from pathlib import Path
from queue import Queue


class QueueHandler(logging.Handler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.log_queue = Queue()

    def emit(self, record: str) -> None:
        self.log_queue.put(record)


class LoggerFactory(object):
    LOG_DIRECTORY = str(os.path.join(Path(__file__).parent.parent.resolve(), "modules", "artifacts", "logs"))
    LOG_FILE_NAME = "log.txt"
    _QUEUE_HANDLER_NAME = "queue_handler"

    @classmethod
    def _create_logger(cls) -> logging.Logger:
        logger = logging.getLogger()
        if logger.hasHandlers():
            # Handler is already initialized
            return logger
        logger.setLevel(logging.INFO)
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        # -------------
        # Stream Handler
        # -------------
        # TODO - un-comment if running only modules.
        logger.addHandler(logging.StreamHandler())      # FIXME
        # -------------
        # Queue Handler
        # -------------
        # FIXME
        # queue_handler = QueueHandler()
        # queue_handler.setFormatter(logging.Formatter(log_format))
        # queue_handler.name = cls._QUEUE_HANDLER_NAME
        # logger.addHandler(queue_handler)
        # -------------
        # File Handler
        # -------------
        # Make sure logs directory exists
        # FIXME
        # Path(cls.LOG_DIRECTORY).mkdir(parents=True, exist_ok=True)
        # # Generate log file
        # log_file_path = os.path.join(cls.LOG_DIRECTORY, cls.LOG_FILE_NAME)
        # log_file_handler = logging.FileHandler(log_file_path, mode="w")
        # log_file_handler.setFormatter(logging.Formatter(log_format))
        # logger.addHandler(log_file_handler)
        return logger

    @staticmethod
    def get_logger() -> logging.Logger:
        """
        A static method called by other modules to initialize logger in
        their own module
        """
        logger = LoggerFactory._create_logger()
        return logger

    @classmethod
    def get_queue_handler(cls, logger: logging.Logger) -> typing.Optional[QueueHandler]:
        for handler in logger.handlers:
            if handler.name == cls._QUEUE_HANDLER_NAME:
                return handler
        return None
