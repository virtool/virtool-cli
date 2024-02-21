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
