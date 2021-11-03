import sys


class IgnoreBrokenPipe(object):
    def __init__(self, stream):
        self.stream = stream

        def ignore_einval(fn):
            from functools import wraps
            @wraps(fn)
            def wrapper(*args, **kwargs):
                try:
                    return fn(*args, **kwargs)
                except OSError as exc:
                    if exc.errno != 22:
                        raise exc
                    else:  # mimicking the default SIGPIPE behavior on Windows
                        sys.exit(1)

            return wrapper

        self.write = ignore_einval(lambda data: self.stream.write(data))
        self.flush = ignore_einval(lambda: self.stream.flush())
