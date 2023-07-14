from pathlib import Path
from rich.console import Console

PACKAGE_DIR = Path(__file__).parents[1]

def run(repo_path: Path):
    """
    Creates an empty repository for a Virtool reference, 
    including workflow callers.

    :param repo_path: Path to a given reference directory
    """
    console = Console()
    caller_path = PACKAGE_DIR / 'callers'
    src_path = repo_path / 'src'
    workflow_path = repo_path / '.github/workflows'
    
    try:
        for path in [src_path, workflow_path]:
            path.mkdir(parents=True, exist_ok=True)
    except:
        console.print(f'Problem creating directories', style="red")
        return
    
    console.print(f'Empty reference repository created at {str(repo_path)}')

    try:
        update_callers(repo_path, caller_path)
    except Exception as e:
        console.print(f'Automatic workflow download failed: {e}')
        return
    
    console.print(f'Caller workflows copied into {str(repo_path)}/.github/workflows/')

def update_callers(repo_path: Path, caller_path: Path):
    """
    Copy caller workflows into repo_path's workflow directory

    :param repo_path: Path to a reference repository
    :param caller_path: Path to a directory containing caller workflows
    """
    if not caller_path.exists():
        raise FileNotFoundError
    
    workflow_path = repo_path / '.github/workflows'
    
    for caller_file in caller_path.glob('*.yml'):
        f = caller_file.read_text()
        caller_path = workflow_path / caller_file.name
        caller_path.write_text(f)