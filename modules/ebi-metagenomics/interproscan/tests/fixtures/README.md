# InterProScan Test Database

## Overview

The InterPro team provides a minimal test database distributed with [InterProScan](https://github.com/ebi-pf-team/interproscan/). This database is used by the test suite to validate core functionality.

## Updating the Database

To download or upgrade the test database, use the helper `download_test_db.sh` script:

```bash
./download_test_db.sh <version> [output_folder]
```

### Parameters

- **`<version>`** (required): The InterProScan version tag (e.g., `5.76-107.0`)
- **`[output_folder]`** (optional): Output directory name. Defaults to `interproscan_db`

### Examples

Download the test database for version 5.76-107.0:

```bash
./download_test_db.sh 5.76-107.0
```

## Output Structure

The script creates the following directory structure:

```
interproscan_db/
└── data/
    └── (test database files)
```
