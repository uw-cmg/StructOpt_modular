class DictionaryObject(dict):
    """A dictionary-like object that allows attribute access (both getting and setting) via an object's dot notation."""
    # Unless you really like python, you probably don't want to worry about this implementation.
    def __init__(self, dictionary):
        self.update(dictionary)
        for k, v in self.items():
            setattr(self, k, DictionaryObject._render(v))

    @staticmethod
    def _render(obj):
        # JSON can only handle a few input types: int/float, str, bool, null, dict, and list
        if isinstance(obj, int) or isinstance(obj, float) or isinstance(obj, str) or isinstance(obj, bool) or obj is None:
            return obj
        elif isinstance(obj, dict):
            return DictionaryObject(obj)
        elif isinstance(obj, list):
            return [DictionaryObject._render(v) for v in obj]

    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value

