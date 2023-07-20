from pathlib import Path

def rm_r(path: Path):
    """
    """
    try:
        for chaff in path.iterdir():
            chaff.unlink()
    except Exception as e:
        return e
    path.rmdir()