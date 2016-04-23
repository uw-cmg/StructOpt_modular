class DictionaryObject(object):
    """A Python representation of a dictionary as an object"""

    def __init__(self, dictionary):
        """
        Args:
            dictionary (dict): a dictionary to be converted into an object -- it cannot have any special types
        """
        for k, v in dictionary.items():
            setattr(self, k, DictionaryObject._render(v))

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


    def to_json(self, **kwargs):
        """Gets a JSON representation of the object
        return    str    a JSON representation of the object
        """
        dictionary = {k: v for k, v in self.__dict__.items()}
        default = kwargs.pop("default", lambda o: {k: v for k, v in o.__dict__.items()})
        sort_keys = kwargs.pop("sort_keys", True)
        indent = kwargs.pop("indent", 4)
        return json.dumps(dictionary, default=default, sort_keys=sort_keys, indent=indent, **kwargs)

    def __str__(self):
        return self.to_json()

    def __repr__(self):
        return "{class_}({dict_})".format(class_=self.__class__.__name__, dict_=self.__dict__)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__ if other else False

    def __ne__(self, other):
        return self.__dict__ != other.__dict__ if other else True

    def __hash__(self):
        return hash(id(self))