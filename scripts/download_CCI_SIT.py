from ftplib import FTP
import os
FTP_HOST   = "ftp.awi.de"
D_FTP_root = "/sea_ice/projects/cci/crdp/v4p0-preview1"
D_local    = "/g/data/gv90/da1339/SeaIce"
def download_ftp_tree(ftp: FTP, remote_dir: str, local_dir: str):
    os.makedirs(local_dir, exist_ok=True)
    ftp.cwd(remote_dir)
    entries = ftp.nlst()
    for entry in entries:
        try:
            ftp.cwd(entry)
            download_ftp_tree(ftp, f"{remote_dir}/{entry}", f"{local_dir}/{entry}")
            ftp.cwd("..")
        except Exception:
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
    ftp.login()
    download_ftp_tree(ftp, D_FTP_root, D_local)
    ftp.quit()
if __name__ == "__main__":
    main()
