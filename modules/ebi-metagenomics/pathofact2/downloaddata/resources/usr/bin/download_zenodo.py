#!/usr/bin/env python3
# Copyright 2024-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import json
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def download_file(url, output_file, user_agent, max_retries=3, retry_delay=5):
    """Download file using urllib with retry logic"""
    for attempt in range(max_retries):
        try:
            req = Request(url, headers={"User-Agent": user_agent})
            with urlopen(req, timeout=3600) as response:
                with open(output_file, "wb") as f:
                    f.write(response.read())
            return True

        except (URLError, HTTPError) as e:
            if attempt < max_retries - 1:
                print(f"Download attempt {attempt + 1} failed: {e}", file=sys.stderr)
                print(f"Retrying in {retry_delay} seconds...", file=sys.stderr)
                time.sleep(retry_delay)
            else:
                print(
                    f"Error downloading file after {max_retries} attempts: {e}",
                    file=sys.stderr,
                )
                return False
        except Exception as e:
            print(f"Unexpected error downloading file: {e}", file=sys.stderr)
            return False

    return False


def get_zenodo_metadata(zenodo_id, user_agent):
    """Fetch Zenodo record metadata using urllib"""
    api_url = f"https://zenodo.org/api/records/{zenodo_id}"

    try:
        req = Request(api_url, headers={"User-Agent": user_agent})
        with urlopen(req, timeout=30) as response:
            data = response.read().decode("utf-8")
            return json.loads(data)
    except HTTPError as e:
        print(
            f"HTTP Error fetching Zenodo metadata: {e.code} - {e.reason}",
            file=sys.stderr,
        )
        sys.exit(1)
    except URLError as e:
        print(f"URL Error fetching Zenodo metadata: {e.reason}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON response: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error fetching metadata: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Download files from Zenodo records")
    parser.add_argument("zenodo_id", type=int, help="Zenodo record ID to download")
    parser.add_argument(
        "--output-file",
        "-o",
        help="Output filename (default: uses original filename from Zenodo)",
    )
    parser.add_argument(
        "--file-index",
        type=int,
        default=0,
        help="File index to download if record has multiple files (default: 0)",
    )
    parser.add_argument(
        "--user-agent",
        default="Python/Zenodo-Downloader",
        help="User agent string for requests",
    )

    args = parser.parse_args()

    zenodo_id = args.zenodo_id
    user_agent = args.user_agent

    # Fetch metadata
    print(f"Fetching Zenodo record metadata for ID: {zenodo_id}")
    api_data = get_zenodo_metadata(zenodo_id, user_agent)

    # Get version info
    db_version = api_data.get("metadata", {}).get("version", "unknown")
    print(f"Record version: {db_version}")

    # Get download URL
    if not api_data.get("files"):
        print("Error: No files found in Zenodo record", file=sys.stderr)
        sys.exit(1)

    # Check if requested file index exists
    if args.file_index >= len(api_data["files"]):
        print(
            f"Error: File index {args.file_index} not found. Record has {len(api_data['files'])} file(s)",
            file=sys.stderr,
        )
        sys.exit(1)

    file_entry = api_data["files"][args.file_index]
    filename = file_entry["key"]
    file_size = file_entry.get("size", "unknown")
    download_url = f"https://zenodo.org/records/{zenodo_id}/files/{filename}"

    # Determine output filename
    output_file = args.output_file if args.output_file else filename

    print(f"Downloading: {filename}")
    print(f"Size: {file_size} bytes" if file_size != "unknown" else "Size: unknown")
    print(f"URL: {download_url}")
    print(f"Saving to: {output_file}")

    # Download file
    if not download_file(download_url, output_file, user_agent):
        print("Error: Failed to download file", file=sys.stderr)
        sys.exit(1)

    print(f"Download complete: {output_file}")


if __name__ == "__main__":
    main()
