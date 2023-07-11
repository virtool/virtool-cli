import json
from pathlib import Path

def get_cache_dir(secret_path: Path):
    with open(secret_path) as f:
        secrets = json.load(f)

    return Path(secrets['repo_dir']) / secrets['cache_dir']

def get_github_token(secret_path: Path):
    with open(secret_path) as f:
        secrets = json.load(f)

    return secrets['GH_BOT_PAT']