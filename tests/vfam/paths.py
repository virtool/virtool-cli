from pathlib import Path

TEST_FILES_PATH = Path(__file__).parents[1] / "files"
VFAM_INPUT_PATH = TEST_FILES_PATH / "vfam_input"
VFAM_INTERMEDIATES_PATH = TEST_FILES_PATH / "vfam_intermediates"

if __name__ == '__main__':
    print(TEST_FILES_PATH)