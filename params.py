class Params(object):
    """
    Object to hold arbitrary attributes with description.

    usage
    -----

    ```
    params = Params('string to be used as title')
    params
    > 'string to be used as title'
    params.any_atttribute = 42
    params.another_config_attr = [777, 'nice']
    ```
    """

    def __init__(self, *args):
        self.__header__ = str(args[0]) if args else None

    def __repr__(self):
        if self.__header__ is None:
             return super(Params, self).__repr__()
        return self.__header__

    def __next__(self):
        """ Fake iteration functionality.
        """
        raise StopIteration

    def __iter__(self):
        """ Fake iteration functionality.
        We skip magic attribues and Structs, and return the rest.
        """
        ks = self.__dict__.keys()
        for k in ks:
            if not k.startswith('__') and not isinstance(k, Params):
                yield getattr(self, k)

    def __len__(self):
        """ Don't count magic attributes or Structs.
        """
        ks = self.__dict__.keys()
        return len([k for k in ks if not k.startswith('__')\
                    and not isinstance(k, Params)])
