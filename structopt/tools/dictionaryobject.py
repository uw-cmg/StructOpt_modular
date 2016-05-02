class DictionaryObject(dict):
    def __init__(self, dictionary):
        self.update(dictionary)
        for key, value in self.items():
            if isinstance(value, dict):
                setattr(self, key, DictionaryObject(value))

    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value

