#!/usr/bin/env python3
"""Convert Jupyter notebooks in src/ to R Markdown files in rmd/."""

import json
import os
import re
import sys

SRC_DIR = os.path.join(os.path.dirname(__file__), "src")
OUT_DIR = os.path.join(os.path.dirname(__file__), "rmd")

COLAB_SKIP_PATTERNS = [
    "from google.colab import",
    "import google.colab",
    "drive.mount(",
    "condacolab",
    "%load_ext rpy2.ipython",
]

R_MAGIC_RE = re.compile(r"^%%R(\s.*)?$", re.MULTILINE)
BASH_MAGIC_RE = re.compile(r"^%%bash\s*$", re.MULTILINE)


def is_colab_cell(source: str) -> bool:
    return any(p in source for p in COLAB_SKIP_PATTERNS)


def prettify_title(stem: str) -> str:
    """Turn a filename stem into a human-readable title."""
    # Replace & and _ with spaces, strip leading numbers
    title = stem.replace("&", "and").replace("_", " ")
    title = re.sub(r"^\d+\s*\d*\s*", "", title).strip()
    return title.title() if title else stem


def safe_stem(name: str) -> str:
    """Make a filesystem-safe stem from notebook name."""
    stem = os.path.splitext(name)[0]
    return stem.replace("&", "_").replace(" ", "_")


def convert_code_cell(source: str) -> str:
    """Convert a code cell source to an Rmd code chunk string."""
    lines = source.splitlines(keepends=True)
    if not lines:
        return ""

    first = lines[0].rstrip("\n").rstrip()

    # %%R or %%R -i var1 -i var2 ...
    if re.match(r"^%%R(\s|$)", first):
        lang = "r"
        body = "".join(lines[1:])
    # %%bash
    elif first == "%%bash":
        lang = "bash"
        body = "".join(lines[1:])
    else:
        lang = "python"
        body = source

    body = body.rstrip()
    if not body.strip():
        return ""

    return f"```{{{lang}}}\n{body}\n```\n"


def notebook_to_rmd(nb_path: str, out_path: str) -> None:
    with open(nb_path, encoding="utf-8") as f:
        nb = json.load(f)

    stem = os.path.splitext(os.path.basename(nb_path))[0]
    title = prettify_title(stem)

    parts = [
        f'---\ntitle: "{title}"\noutput: html_document\n---\n\n'
    ]

    for cell in nb.get("cells", []):
        cell_type = cell.get("cell_type", "")
        source = "".join(cell.get("source", []))

        if not source.strip():
            continue

        if cell_type == "markdown":
            parts.append(source.rstrip() + "\n\n")
        elif cell_type == "code":
            if is_colab_cell(source):
                continue
            chunk = convert_code_cell(source)
            if chunk:
                parts.append(chunk + "\n")

    rmd_content = "".join(parts)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(rmd_content)

    print(f"  {os.path.basename(nb_path)} -> {os.path.relpath(out_path)}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    notebooks = sorted(
        f for f in os.listdir(SRC_DIR) if f.endswith(".ipynb")
    )
    if not notebooks:
        print("No notebooks found in src/", file=sys.stderr)
        sys.exit(1)

    print(f"Converting {len(notebooks)} notebooks to {OUT_DIR}/")
    for nb_name in notebooks:
        nb_path = os.path.join(SRC_DIR, nb_name)
        out_name = safe_stem(nb_name) + ".Rmd"
        out_path = os.path.join(OUT_DIR, out_name)
        notebook_to_rmd(nb_path, out_path)

    print("Done.")


if __name__ == "__main__":
    main()
