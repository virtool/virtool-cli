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

