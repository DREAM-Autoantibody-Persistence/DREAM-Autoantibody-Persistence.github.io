"""Download the eight documented public DREAM workbooks, verifying SHA256."""
import argparse
import csv
import hashlib
import ssl
import time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import URLError

PROJECT_DIR = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=PROJECT_DIR / "data")
    args = parser.parse_args()
    args.data_dir.mkdir(parents=True, exist_ok=True)
    with (PROJECT_DIR / "source_manifest.csv").open(encoding="utf-8", newline="") as handle:
        sources = list(csv.DictReader(handle))
    context = ssl.create_default_context()
    for item in sources:
        target = args.data_dir / item["File"]
        if target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == item["Hash"]:
            print("Already verified " + target.name)
            continue
        temporary = target.with_suffix(".download")
        try:
            for attempt in range(3):
                try:
                    request = Request(item["SourceUrl"], headers={"User-Agent": "DREAM-reproducible-analysis/1.0"})
                    with urlopen(request, timeout=120, context=context) as response, temporary.open("wb") as output:
                        while True:
                            chunk = response.read(1024*1024)
                            if not chunk:
                                break
                            output.write(chunk)
                    break
                except (URLError, TimeoutError, ConnectionError):
                    if attempt == 2:
                        raise
                    print("Retrying " + target.name)
                    time.sleep(attempt+1)
            if hashlib.sha256(temporary.read_bytes()).hexdigest() != item["Hash"]:
                raise ValueError("Downloaded file differs from the frozen source: " + target.name)
            temporary.replace(target)
            print("Downloaded and verified " + target.name)
        finally:
            if temporary.exists():
                temporary.unlink()

if __name__ == "__main__":
    main()
