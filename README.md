# virtool-cli

A command line tool for working with Virtool data.

## Installation

```shell script
pip install virtool-cli
```

## Usage

### Ref
Commands related to building, maintaining and pulling new data for reference databases.

#### Build

To build a `reference.json` file at **OUTPUT_PATH** from a src directory at **DIRECTORY_PATH**:

```shell script
virtool ref build -src DIRECTORY_PATH -o OUTPUT_PATH
```

For a more easily readable output file, append the `--indent` flag:

```shell script
virtool ref build -src DIRECTORY_PATH -o OUTPUT_PATH --indent
```

To specify a version to include in the `reference.json` file:

```shell script
virtool ref build -src DIRECTORY_PATH --version VERSION
```

#### Divide

To split a single JSON file reference at **INPUT_PATH**
into a readable hierarchal directory at **OUTPUT_PATH**:

```shell script
virtool ref divide --file_path INPUT_PATH -o OUTPUT_PATH
```

#### Command subgroup: Migrate

To update an existing reference source directory at **DIRECTORY_PATH** generated by virtool-cli v1:

```shell script
virtool ref migrate -src DIRECTORY_PATH
```

#### Command subgroup: Update

##### Fetch updates for a single OTU

```shell script
virtool ref update otu -otu OTU_PATH --catalog_path CATALOG_PATH
```

##### Fetch updates for all OTUs in a reference

```shell script
virtool ref update reference -src SOURCE_PATH --catalog_path CATALOG_PATH
```

### Environmental Variables
Some of the tools in the CLI make API requests to NCBI. 
Unauthenticated requests are are limited to 3 per second. 
Setting NCBI account credentials in your environmental variables can 
increase the allowed rate to 10 requests per second.

| Name | Description |
|----|---------|
| `NCBI_EMAIL` | The e-mail address used for your NCBI account |
| `NCBI_API_KEY` | The [API key](https://www.ncbi.nlm.nih.gov/account/settings/) associated with your NCBI account. |

## Tests

We prefer running the test suite in Docker.

1. Make sure the image builds:
   ```
   docker build -q --target test .
   ```
   
2. Run Pytest in the container:
   ```
   docker run --rm -t $(docker build -q --target test .) pytest
   ```
