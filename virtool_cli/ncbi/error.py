class NCBIClientError(Exception):
    def __init__(self, message: str = ""):
        self.message = message


class NCBIParseError(Exception):
    def __init__(self, keys, message: str = ""):
        self.keys = list(keys)
        self.message = message
