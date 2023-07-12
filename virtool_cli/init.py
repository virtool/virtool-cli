import requests
import tarfile
from pathlib import Path
from typing import Optional
import structlog


logger = structlog.get_logger()

def run(repo_path: Path, api_url: str, cache_path: Path, gh_token=''):
    """
    Creates an empty repository for a Virtool reference, 
    including workflow callers.
    Workflow callers are downloaded from the latest release URL 
    of the workflow repository.

    :param ref_dir: Path to a given reference directory
    :param api_url: GitHub API URL of the latest caller.tar.gz release asset
    :param cache_path: Path to a cache directory
    :param gh_token: A GitHub token for use with private repositories (Optional)
    """
    src_path = repo_path / 'src'
    workflow_path = repo_path / '.github/workflows'
    
    try:
        for path in [src_path, workflow_path]:
            path.mkdir(parents=True, exist_ok=True)
    except:
        logger.error('Problem creating directories')
        return
    
    logger.info(
        'Empty reference repository created', 
        repo_path=str(repo_path))

    try:
        caller_tar = fetch_callers(
            api_url=api_url,
            cache_path=cache_path, 
            gh_token=gh_token)
    except Exception as e:
        logger.error('%s', e, api_url=api_url)
        return

    try:
        update_callers(caller_tar, cache_path, repo_path)
    except Exception as e:
        logger.error('Automatic workflow download failed: %s', e)
        return
    
    logger.debug(
        'Caller workflows downloaded into .github/workflows/', 
        repo_path=str(repo_path))

def fetch_callers(api_url: str, cache_path: Path, gh_token='') -> Optional[Path]:
    """
    Download caller library to cache

    :param api_url: GitHub API URL of the latest caller.tar.gz release asset
    :param cache_path: Path to a cache directory
    :param gh_token: A GitHub token for use with private repositories (Optional)
    """
    headers = { 'Accept': 'application/vnd.github+json' }
    
    if gh_token:
        headers['Authorization'] = f'token {gh_token}'

    r = requests.get(api_url, headers=headers)
    if r.status_code != 200:
        raise requests.HTTPError(
            f'HTTP Error {r.status_code}: Could not find file to download at api_url')

    try:
        r_assets = requests.get(r.json()['assets_url'], headers=headers)
        caller_url = r_assets.json()[0]['url']
    except (requests.JSONDecodeError, KeyError):
        raise KeyError(
            f'Error: Could not extract a download url from api_url')
    
    headers['Accept'] = 'application/octet-stream'
    r_callers = requests.get(caller_url, headers=headers)

    try:
        tar_path = cache_path / 'callers.tar.gz'
        with open(tar_path, "wb") as f:
            for chunk in r_callers.iter_content(chunk_size=512):
                if chunk: f.write(chunk)
    except:
        raise requests.HTTPError(f'File download to {tar_path} failed')

    if not tar_path.exists():
        raise FileNotFoundError(f'File download to {tar_path} failed')
    
    try:
        caller_data = tarfile.open(tar_path, mode="r:gz")
        caller_data.extractall(path=cache_path)
        caller_data.close()

        cached_caller_path = cache_path / 'callers'
        if not cached_caller_path.exists():
            raise tarfile.ExtractError(f'Encountered a problem with tarfile at {tar_path}')
    except:
        raise tarfile.ExtractError(f'Encountered a problem with tarfile at {tar_path}')

    return tar_path

def update_callers(tar_path: Path, cache_path: Path, repo_path: Path):
    """
    Expand downloaded caller archive into cache
    and copy contents into repo_path's workflow directory

    :param tar_path: Path to a cached archive of caller workflows 
    :param cache_path: Path to a cache directory
    :param repo_path: Path to a reference repository
    """

    cached_caller_path = cache_path / 'callers'
    if not cached_caller_path.exists():
        raise tarfile.ExtractError(f'Encountered a problem with tarfile at {tar_path}')
    
    workflow_path = repo_path / '.github/workflows'
    
    for caller_file in cached_caller_path.glob('*.yml'):
        f = caller_file.read_text()
        caller_path = workflow_path / caller_file.name
        caller_path.write_text(f)