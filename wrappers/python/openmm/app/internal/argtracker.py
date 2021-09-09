from inspect import signature

class ArgTracker(dict):
    """This class tracks usage of function keyword arguments to detect ones that are never used."""

    def __init__(self, *args):
        super().__init__(*args)
        self.accessed = set()

    def __getitem__(self, key):
        self.accessed.add(key)
        return super().__getitem__(key)

    def get(self, key, default):
        self.accessed.add(key)
        return super().get(key, default)

    def checkArgs(self, fn):
        """Throw an exception if any argument was never used."""
        parameters = signature(fn).parameters
        for key in self:
            if key not in self.accessed:
                if key not in parameters or self[key] != parameters[key].default:
                    raise ValueError(f"The argument '{key}' was specified to {fn.__name__}() but was never used.")
