class FrozenAttributesTrait:
    """
    Class which prevents adding new attributes after a call to freeze.
    This should prevent us from doing (some?) typos in our config files.
    Based on this stackoverflow post https://stackoverflow.com/a/29368642.
    """
    __attributes_frozen = False
    def __setattr__(self, key, value):
        if self.__attributes_frozen and not hasattr(self, key):
            raise TypeError( "attribute {} cannot be added to {} ".format(key, self) )
        object.__setattr__(self, key, value)

    def _freeze_attributes(self):
        self.__attributes_frozen = True
