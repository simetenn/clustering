import subprocess
import os
import argparse


def convert(filename, output_dir):
    cmd = ["libreoffice", "--headless", "--convert-to", "csv"] + [filename] + ["--outdir"] + [output_dir]
    subprocess.Popen(cmd)


def convertFolder(foldername):
    for root, dirs, files in os.walk(foldername):
        for filename in files:
            if filename.endswith(".xls"):
                convert(os.path.join(root, filename), root)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert all xls files to csv")
    parser.add_argument("foldername")

    args = parser.parse_args()

    convertFolder(args.foldername)
