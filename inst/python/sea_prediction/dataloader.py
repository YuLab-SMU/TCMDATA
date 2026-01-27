import os
import gdown

def download_data(file_path, file_id):
    """
    Download database and parameters files from Google Drive.
    Checks if a file exists locally. If not, downloads it from Google Drive using its File ID.
    """
    if not os.path.exists(file_path):
        print(f"Downloading core data file: {file_path} ...")
        # Construct the direct download URL for Google Drive
        url = f'https://drive.google.com/uc?id={file_id}'
        try:
            # gdown handles large file warnings automatically
            gdown.download(url, file_path, quiet=False)
        except Exception as e:
            print(f"Failed to download {file_path}: {e}")
    else:
        print(f"File already exists: {file_path}")

if __name__ == "__main__":
    # --- Configuration: Users can modify these parameters ---
    
    # Google Drive File IDs (Extracted from shareable links)
    DB_ID = "1SV7-vrG3BFL51Ie07jtEBaOlFv8pkK3D"
    FIT_ID = "1jviqHWp0UaXj2bvhO5WY_bnOUyqLs_mu"
    
    # Local destination paths
    DB_PATH = "target_fp_db.pkl"
    FIT_PATH = "chembl27_rdkit_ecfp4.fit"

    # Run download checks
    download_data(DB_PATH, DB_ID)
    download_data(FIT_PATH, FIT_ID)