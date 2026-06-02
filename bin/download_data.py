#!/usr/bin/env python
"""
Download large data files for streamobs.

This script downloads a data archive (zip file) from a remote repository
(e.g., Zenodo, institutional server) and extracts it to the data/ directory.
"""

import os
import sys
import argparse
import urllib.request
import urllib.error
from pathlib import Path
import zipfile
import tempfile


# =============================================================================
# CONFIGURATION - Update this URL when data location changes
# =============================================================================

# Base URL where data archive is hosted
BASE_DATA_URL = "https://zenodo.org/records/18298544/files/"

# Name of the data archive file (should be a zip file)
DATA_ARCHIVE_NAME = "data.zip"

# Full URL to download
DATA_ARCHIVE_URL = BASE_DATA_URL + DATA_ARCHIVE_NAME

# Expected size (approximate, for user information)
ARCHIVE_SIZE_MB = 30 #Mb

# =============================================================================
# =============================================================================


def download_file(url, output_path, description="file"):
    """Download a file with progress indication."""
    print(f"Downloading {description}...")
    print(f"  From: {url}")
    print(f"  To: {output_path}")
    
    def progress_hook(count, block_size, total_size):
        """Show download progress."""
        if total_size > 0:
            percent = int(count * block_size * 100 / total_size)
            mb_downloaded = count * block_size / (1024 * 1024)
            mb_total = total_size / (1024 * 1024)
            sys.stdout.write(f"\r  Progress: {percent}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)")
            sys.stdout.flush()
    
    try:
        urllib.request.urlretrieve(url, output_path, reporthook=progress_hook)
        print("\n  ✓ Download complete!")
        return True
    except Exception as e:
        print(f"\n  ✗ Download failed: {e}")
        return False


def extract_zip(zip_path, extract_to, description="archive"):
    """Extract a zip file."""
    print(f"\nExtracting {description}...")
    print(f"  From: {zip_path}")
    print(f"  To: {extract_to}")
    
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Get list of files
            file_list = zip_ref.namelist()
            print(f"  Found {len(file_list)} files in archive")
            
            # Extract all files
            zip_ref.extractall(extract_to)
            print("  ✓ Extraction complete!")
            
        # Clean up unwanted files after extraction
        cleanup_unwanted_files(extract_to)
        
        return True
    except Exception as e:
        print(f"  ✗ Extraction failed: {e}")
        return False


def cleanup_unwanted_files(base_path):
    """Remove unwanted files like .DS_Store, .backup, etc. after extraction."""
    base_path = Path(base_path)
    removed_count = 0
    
    # Patterns to remove
    unwanted_patterns = [
        '**/.DS_Store',      # macOS system files
        '**/__MACOSX',       # macOS resource forks
        '**/*.backup',       # Backup files
        '**/*.bak',          # Backup files
        '**/*~',             # Temporary files
    ]
    
    for pattern in unwanted_patterns:
        for item in base_path.glob(pattern):
            try:
                if item.is_file():
                    item.unlink()
                    removed_count += 1
                elif item.is_dir():
                    import shutil
                    shutil.rmtree(item)
                    removed_count += 1
            except Exception:
                pass  # Ignore errors during cleanup
    
    if removed_count > 0:
        print(f"  Cleaned up {removed_count} unwanted file(s)")



def list_data_contents(data_dir):
    """List what's in the data directory after download."""
    print("\n📂 Data directory contents:")
    print("=" * 80)
    
    if not data_dir.exists():
        print("  (empty - data directory doesn't exist yet)")
        return
    
    # List subdirectories
    subdirs = [d for d in data_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
    if not subdirs:
        print("  (empty - no subdirectories found)")
        return
    
    for subdir in sorted(subdirs):
        print(f"\n  {subdir.name}/")
        # Count files in subdirectory
        try:
            files = list(subdir.rglob('*'))
            files = [f for f in files if f.is_file() and not f.name.startswith('.')]
            total_size = sum(f.stat().st_size for f in files)
            print(f"    Files: {len(files)} ({total_size / (1024*1024):.1f} MB)")
        except Exception:
            print(f"    (unable to read directory)")



def main():
    parser = argparse.ArgumentParser(
        description='Download and extract large data files for streamobs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download and extract all data
  python download_data.py
  
  # Force re-download even if data exists
  python download_data.py --force
  
  # List what's currently in the data directory
  python download_data.py --list
  
  # Use custom data URL
  python download_data.py --url https://my-server.edu/data.zip
        """
    )
    
    parser.add_argument('--list', action='store_true',
                        help='List current data directory contents without downloading')
    parser.add_argument('--url', type=str, default=DATA_ARCHIVE_URL,
                        help=f'URL for data archive (default: configured Zenodo URL)')
    parser.add_argument('--data-dir', type=str, default=None,
                        help='Data directory (default: streamobs/data/)')
    parser.add_argument('--force', action='store_true',
                        help='Force re-download even if data exists')
    parser.add_argument('--keep-archive', action='store_true',
                        help='Keep the downloaded zip file after extraction')
    
    args = parser.parse_args()
    
    # Determine data directory
    if args.data_dir:
        data_dir = Path(args.data_dir)
    else:
        # Assume script is in streamobs/bin/
        script_dir = Path(__file__).parent
        data_dir = script_dir.parent / 'data'
    
    # List contents and exit if requested
    if args.list:
        list_data_contents(data_dir)
        return 0
    
    # Check if URL is configured
    if 'XXXXX' in args.url:
        print("=" * 80)
        print("ERROR: Data URL is not yet configured!")
        print("=" * 80)
        print("\nThe DATA_ARCHIVE_URL in this script is still set to a placeholder.")
        print("Please update it with the actual data hosting location.\n")
        print("Steps:")
        print("  1. Upload your data/ directory as a zip file to Zenodo or another host")
        print("  2. Edit this script and update BASE_DATA_URL and DATA_ARCHIVE_NAME")
        print("  3. Or use --url to specify a custom URL\n")
        return 1
    
    # Check if data already exists
    if data_dir.exists() and not args.force:
        subdirs = [d for d in data_dir.iterdir() if d.is_dir() and not d.name.startswith('.')]
        if subdirs:
            print("=" * 80)
            print("Data directory already exists with content!")
            print("=" * 80)
            list_data_contents(data_dir)
            print("\n" + "=" * 80)
            print("Use --force to re-download and overwrite, or --list to view contents")
            return 0
    
    # Download and extract
    print("=" * 80)
    print("Stream Simulation Data Download")
    print("=" * 80)
    print(f"\nData URL: {args.url}")
    print(f"Destination: {data_dir}")
    print(f"Archive size: ~{ARCHIVE_SIZE_MB} MB")
    print("\n" + "=" * 80)
    
    # Create data directory
    data_dir.mkdir(parents=True, exist_ok=True)
    
    # Create temporary file for download
    with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as tmp_file:
        tmp_path = Path(tmp_file.name)
    
    try:
        # Download
        if not download_file(args.url, tmp_path, "data archive"):
            print("\n✗ Download failed!")
            return 1
        
        # Extract
        if not extract_zip(tmp_path, data_dir.parent, "data archive"):
            print("\n✗ Extraction failed!")
            return 1
        
        # Show what was extracted
        list_data_contents(data_dir)
        
        print("\n" + "=" * 80)
        print("✓ Data download and extraction complete!")
        print("=" * 80)
        print("\nYou can now run stream simulations with this data.")
        print("The data has been extracted to:", data_dir)
        
        return 0
        
    finally:
        # Clean up temporary file unless requested to keep it
        if tmp_path.exists():
            if args.keep_archive:
                archive_path = data_dir.parent / DATA_ARCHIVE_NAME
                tmp_path.rename(archive_path)
                print(f"\nArchive saved to: {archive_path}")
            else:
                tmp_path.unlink()


if __name__ == '__main__':
    sys.exit(main())
