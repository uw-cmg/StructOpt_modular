import json


class DictionaryObject(dict):
    def __init__(self, dictionary):
        self.update(dictionary)
        for k, v in self.items():
            setattr(self, k, DictionaryObject._render(v))

    @staticmethod
    def _render(obj):
        if isinstance(obj, int) or isinstance(obj, float):
            return obj
        elif  isinstance(obj, str):
            try:
                return eval(obj)
            except:
                return obj
        elif isinstance(obj, dict):
            return DictionaryObject(obj)
        elif isinstance(obj, list):
            return [DictionaryObject._render(v) for v in obj]

    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        self[key] = value

    def to_json(self, *args, **kwargs):
        return json.dumps(self, *args, **kwargs)

    def __str__(self):
        return self.to_json()
