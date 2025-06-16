from ftplib import FTP
import os

# FTP configuration
FTP_HOST = "ftp.awi.de"
FTP_ROOT_DIR = "/sea_ice/projects/cci/crdp/v4p0-preview1"
LOCAL_SAVE_DIR = "/g/data/gv90/da1339/SeaIce"

def download_ftp_tree(ftp: FTP, remote_dir: str, local_dir: str):
    """Recursively download all .nc files from remote_dir into local_dir"""
    os.makedirs(local_dir, exist_ok=True)
    ftp.cwd(remote_dir)
    entries = ftp.nlst()

    for entry in entries:
        try:
            ftp.cwd(entry)  # if this succeeds, it's a directory
            download_ftp_tree(ftp, f"{remote_dir}/{entry}", f"{local_dir}/{entry}")
            ftp.cwd("..")
        except Exception:
            # assume it's a file
            if entry.endswith(".nc"):
                local_file_path = os.path.join(local_dir, entry)
                if not os.path.exists(local_file_path):
                    print(f"Downloading: {remote_dir}/{entry}")
                    with open(local_file_path, "wb") as f:
                        ftp.retrbinary(f"RETR {entry}", f.write)
                else:
                    print(f"Already exists: {local_file_path}")

def main():
    ftp = FTP(FTP_HOST)
    ftp.login()  # anonymous login
    download_ftp_tree(ftp, FTP_ROOT_DIR, LOCAL_SAVE_DIR)
    ftp.quit()

if __name__ == "__main__":
    main()
