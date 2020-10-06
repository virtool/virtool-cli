# virtool-cli

A command line tool for working with Virtool data.

[![Build Status](https://cloud.drone.io/api/badges/virtool/virtool-cli/status.svg)](https://cloud.drone.io/virtool/virtool-cli)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/f4d6416f3d434c62af89c2ba94f8343a)](https://www.codacy.com/gh/virtool/virtool-cli/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=virtool/virtool-cli&amp;utm_campaign=Badge_Grade)

## Installation

```shell script
git clone https://github.com/virtool/virtool-cli.git
pip install .
```

### Usage
#### Build
To build a `reference.json` file from a src directory

```shell script
virtool build -src DIRECTORY_PATH -o OUTPUT_PATH
```

If you wish for the output file to be more easily readable you can specify it to be indented

```shell script
virtool build -src DIRECTORY_PATH -i
```

To specify a version to include in the `reference.json` file

```shell script
virtool build -src DIRECTORY_PATH -V VERSION
```

#### Divide
If you have an existing `reference.json` you can divide it into a src directory

```shell script
virtool divide -src REFERENCE_PATH -o OUTPUT_PATH
```

#### Taxid
To fetch taxon ids for OTUs in a given src directory that do not currently have one

```shell script
virtool taxid -src DIRECTORY_PATH
```

If you wish to force a taxon id lookup for all OTUs

```shell script
virtool taxid -src DIRECTORY_PATH -f
```

### Environmental Variables
Some of the tools in the CLI make requests to NCBI databases, which are usually
rate-limited to 3 requests/second. If you are working with larger directories this will cause some tools
to run extremely slowly. Fortunately, with correct authenication this limit can be raised to
10 requests/second. This involves setting up the following environmental variables:

`NCBI_EMAIL` - The e-mail address used for your NCBI account

`NCBI_API_KEY` - The API key associated with your NCBI account. If you are currently signed into NCBI it can be found [here](https://www.ncbi.nlm.nih.gov/account/settings/).