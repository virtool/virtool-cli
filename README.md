# virtool-cli

A command line tool for working with Virtool data.

[![Build Status](https://cloud.drone.io/api/badges/virtool/virtool-cli/status.svg)](https://cloud.drone.io/virtool/virtool-cli)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f4d6416f3d434c62af89c2ba94f8343a)](https://www.codacy.com/gh/virtool/virtool-cli/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=virtool/virtool-cli&amp;utm_campaign=Badge_Grade)

## Installation

```shell script
pip install virtool-cli
```

## Usage

### Ref
Commands related to references

#### Build
To build a `reference.json` file from a src directory

```shell script
virtool ref build -src DIRECTORY_PATH -o OUTPUT_PATH
```

If you wish for the output file to be more easily readable you can specify it to be indented

```shell script
virtool ref build -src DIRECTORY_PATH -i
```

To specify a version to include in the `reference.json` file

```shell script
virtool ref build -src DIRECTORY_PATH -V VERSION
```

#### Divide
If you have an existing `reference.json` you can divide it into a src directory

```shell script
virtool ref divide -src REFERENCE_PATH -o OUTPUT_PATH
```

#### Taxid
To fetch taxon ids for OTUs in a given src directory that do not currently have one

```shell script
virtool ref taxid -src DIRECTORY_PATH
```

If you wish to force a taxon id lookup for all OTUs

```shell script
virtool ref taxid -src DIRECTORY_PATH -f
```

#### Merge-Refs
To interactively merge two references

```shell script
virtool ref merge -s SOURCE_DIRECTORY_PATH -t TARGET_DIRECTORY_PATH
```

If you wish to resume an incomplete merge

```shell script
virtool ref merge -s SOURCE_DIRECTORY_PATH -t TARGET_DIRECTORY_PATH -r
```

If you wish to place new isolates directly into target reference

```shell script
virtool ref merge -s SOURCE_DIRECTORY_PATH -t TARGET_DIRECTORY_PATH -i
```

If you wish to specify a name for the merged reference directory

```shell script
virtool ref merge -s SOURCE_DIRECTORY_PATH -t TARGET_DIRECTORY_PATH -o OUTPUT_DIRECTORY_NAME
```

If you wish to produce a CSV file containing information about all isolates with unknown source name or source type

```shell script
virtool ref merge -s SOURCE_DIRECTORY_PATH -t TARGET_DIRECTORY_PATH --output-unknown-isolates
```

If you wish to specify the name of the CSV file containing all isolates with unknown source name or source type

```shell script
virtool ref merge -s SOURCE_DIRECTORY_PATH -t TARGET_DIRECTORY_PATH --output-unknown-isolates --unknown-output OUTPUT_FILE_NAME
```

### Environmental Variables
Some of the tools in the CLI make API requests to NCBI. Unauthenticated requests are are limited to 3 per second. Setting NCBI credentials in environmental variables can increase this to 10 per second.

| Name | Description |
|----|---------|
| `NCBI_EMAIL` | The e-mail address used for your NCBI account |
| `NCBI_API_KEY` | The [API key](https://www.ncbi.nlm.nih.gov/account/settings/) associated with your NCBI account. |

