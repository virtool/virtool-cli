class NCBIClientError(Exception):
    def __init__(self, message: str = ""):
        self.message = message


class IncompleteRecordsError(NCBIClientError):
    def __init__(self, message: str = "", data: list[dict] = None):
        super().__init__(message)
        if data is None:
            self.data = []
        else:
            self.data = data


class NCBIParseError(Exception):
    def __init__(self, keys, message: str = ""):
        self.keys = list(keys)
        self.message = message
