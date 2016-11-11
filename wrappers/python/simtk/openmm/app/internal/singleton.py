"""
Creates a subclass for all classes intended to be a singleton. This
maintains the correctness of instance is instance even following
pickling/unpickling
"""
class Singleton(object):
    _inst = None
    def __new__(cls):
        if cls._inst is None:
            cls._inst = super(Singleton, cls).__new__(cls)
        return cls._inst

    def __reduce__(self):
        return repr(self)

