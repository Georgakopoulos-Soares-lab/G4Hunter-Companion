if __name__ == "__main__":
    import sys
    from pathlib import Path

    # Ensure the script is run from the correct directory
    suffix = sys.argv[2] if len(sys.argv) > 2 else "fna"
    # Add the parent directory to the system path
    files = [file for file in Path(sys.argv[1]).glob(f"*.{suffix}") if file.is_file()]
    with open("design.csv", mode="w", encoding="utf-8") as f:
        f.write("accession_id\taccession_path\n")
        for file in files:
            f.write(f"{file.stem}\t{file.as_posix()}\n")


