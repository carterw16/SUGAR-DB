# Downloading mdbtools
import os
import platform
import tarfile
from urllib.request import urlretrieve

current_dir = os.path.realpath(os.path.dirname(__file__))
tar_file_name = os.path.join(current_dir, "mdbtools.tar.gz")
mdb_dir = os.path.join(current_dir, "mdbtools")

if platform.system() == "Windows":
    URL = "https://github.com/kdheepak/mdbtools/releases/download/download/mdbtools-windows.tar.gz"
elif platform.system() == "Darwin":
    URL = "https://github.com/kdheepak/mdbtools/releases/download/download/mdbtools-osx.tar.gz"
else:
    URL = "https://github.com/kdheepak/mdbtools/releases/download/download/mdbtools-linux.tar.gz"


def download_mdbtools(url):
    urlretrieve(url, tar_file_name)
    with tarfile.open(tar_file_name) as tf:
        tf.extractall(mdb_dir)
