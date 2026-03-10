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

# Matches the "Open in Colab" badge line in markdown cells
COLAB_BADGE_RE = re.compile(r"colab\.research\.google\.com", re.IGNORECASE)

R_MAGIC_RE = re.compile(r"^%%R(\s.*)?$", re.MULTILINE)
BASH_MAGIC_RE = re.compile(r"^%%(bash|shell)\s*$", re.MULTILINE)
OS_ENVIRON_RE = re.compile(r"os\.environ\[(['\"])(\w+)\1\]\s*=\s*(.+)")


def is_colab_cell(source: str) -> bool:
    return any(p in source for p in COLAB_SKIP_PATTERNS)


def has_os_environ(source: str) -> bool:
    return bool(OS_ENVIRON_RE.search(source))


def python_val_to_bash(expr: str, var_values: dict) -> str:
    """Convert a simple Python expression to a bash-safe quoted value."""
    expr = expr.strip()

    # String literal (single or double quoted)
    if len(expr) >= 2 and expr[0] in ('"', "'") and expr[-1] == expr[0]:
        inner = expr[1:-1]
        return f'"{inner}"'

    # os.getcwd()
    if expr == "os.getcwd()":
        return '"$(pwd)"'

    # Simple variable reference
    if re.match(r'^\w+$', expr):
        if expr in var_values:
            return var_values[expr]
        return f'"${{{expr}}}"'

    # String concatenation (handles var + "str", "str" + var, os.getcwd() + "str", etc.)
    parts = re.split(r'\s*\+\s*', expr)
    if len(parts) > 1:
        bash_inner = ""
        for part in parts:
            part = part.strip()
            if part == "os.getcwd()":
                bash_inner += "$(pwd)"
            elif re.match(r'^\w+$', part):
                if part in var_values:
                    val = var_values[part]
                    # Unquote the stored value for embedding
                    if len(val) >= 2 and val[0] in ('"', "'"):
                        bash_inner += val[1:-1]
                    else:
                        bash_inner += val
                else:
                    bash_inner += f'${{{part}}}'
            elif len(part) >= 2 and part[0] in ('"', "'") and part[-1] == part[0]:
                bash_inner += part[1:-1]
            else:
                bash_inner += part
        return f'"{bash_inner}"'

    # Fallback: wrap in quotes
    return f'"{expr}"'


def convert_os_environ_cell(source: str) -> str:
    """Convert a Python cell with os.environ calls to a bash export chunk."""
    lines = source.splitlines()

    # Collect simple variable assignments: name -> bash-quoted value
    var_values = {}
    for line in lines:
        stripped = line.strip()
        m = re.match(r'^(\w+)\s*=\s*(.+)$', stripped)
        if m and not stripped.startswith('#') and not stripped.startswith('import'):
            varname = m.group(1)
            rhs = m.group(2).strip()
            # Store the bash-equivalent value
            var_values[varname] = python_val_to_bash(rhs, var_values)

    exports = []
    for line in lines:
        m = OS_ENVIRON_RE.search(line.strip())
        if m:
            key = m.group(2)
            rhs = m.group(3).strip()
            value = var_values.get(rhs, python_val_to_bash(rhs, var_values))
            exports.append(f'export {key}={value}')

    if not exports:
        return ""

    body = "\n".join(exports)
    return f"```{{bash}}\n{body}\n```\n"


def prettify_title(stem: str) -> str:
    """Turn a filename stem into a human-readable title."""
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
    # %%bash or %%shell
    elif re.match(r"^%%(bash|shell)\s*$", first):
        lang = "bash"
        body = "".join(lines[1:])
    else:
        lang = "python"
        body = source

    body = body.rstrip()
    if not body.strip():
        return ""

    # Convert os.environ Python cells to bash export chunks
    if lang == "python" and has_os_environ(body):
        return convert_os_environ_cell(body)

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
            # Strip the "Open in Colab" badge line(s)
            filtered_lines = [
                l for l in source.splitlines(keepends=True)
                if not COLAB_BADGE_RE.search(l)
            ]
            source = "".join(filtered_lines).strip()
            if source:
                parts.append(source + "\n\n")
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
