from pathlib import Path

TEST_FILES_PATH = Path(__file__).parent / "files"
VFAM_INPUT_PATH = TEST_FILES_PATH / "input"
VFAM_INTERMEDIATES_PATH = TEST_FILES_PATH / "intermediates"

if __name__ == '__main__':
    print(TEST_FILES_PATH)