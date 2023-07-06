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

#### Repair
Fix folder-JSON name mismatches and incorrect taxid types

```shell script
virtool ref repair -src DIRECTORY_PATH
```

#### taxid
Search GenBank for matching OTUs and add their taxon ids to otu.json entries in the source directory.

```shell script
virtool ref taxid -src DIRECTORY_PATH
```

### Environmental Variables
Some of the tools in the CLI make API requests to NCBI. Unauthenticated requests are are limited to 3 per second. Setting NCBI credentials in environmental variables can increase this to 10 per second.

| Name | Description |
|----|---------|
| `NCBI_EMAIL` | The e-mail address used for your NCBI account |
| `NCBI_API_KEY` | The [API key](https://www.ncbi.nlm.nih.gov/account/settings/) associated with your NCBI account. |

