import time


class Timer(object):
    def __enter__(self):
        self.start_time = time.monotonic()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.elapsed_time = time.monotonic() - self.start_time
