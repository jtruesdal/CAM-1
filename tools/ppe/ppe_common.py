"""Shared helpers for the SCAM PPE workflow (tools/ppe).

Design: one base SCAM case is built and run once with CIME. Each PPE member
is a plain directory holding copies of the base run directory's namelist and
config files, with parameter values patched into atm_in. All members run the
same cesm.exe from their own directory. No CIME clone per member.
"""
import csv
import re
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

# Files copied from the base RUNDIR into every member directory.
RUN_FILES = [
    "atm_in", "drv_in", "drv_flds_in", "lnd_in", "ice_in", "docn_in",
    "docn.streams.xml", "fd.yaml", "nuopc.runconfig", "nuopc.runseq",
]
# Directories created (empty) in every member directory.
RUN_DIRS = ["timing", "timing/checkpoints", "init_generated_files"]


def load_config(path):
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    for key in ("base_rundir", "exe", "env_file", "ppe_root"):
        if key not in cfg:
            sys.exit(f"config: missing required key '{key}'")
    cfg["base_rundir"] = Path(cfg["base_rundir"])
    cfg["ppe_root"] = Path(cfg["ppe_root"])
    return cfg


def member_name(i):
    return f"m{i:04d}"


# ---------------------------------------------------------------- atm_in
# atm_in written by build-namelist has one "name = value" per line, which
# lets us patch values in place and leave every other byte unchanged.
def _var_re(name):
    return re.compile(rf"^(\s*{re.escape(name)}\s*=\s*)(\S[^\n]*?)(\s*)$",
                      re.IGNORECASE | re.MULTILINE)


def read_nml_value(text, name):
    """Return the raw value string of a scalar namelist variable, or None."""
    matches = _var_re(name).findall(text)
    if len(matches) > 1:
        sys.exit(f"atm_in: '{name}' appears {len(matches)} times")
    return matches[0][1] if matches else None


def fortran_to_float(s):
    return float(s.strip().lower().replace("d", "e"))


def format_value(v):
    """Full-precision Fortran real literal."""
    return f"{v:.15e}".replace("e", "D")


def patch_nml(text, values):
    """Replace values of existing scalar variables. Refuses to add new ones."""
    for name, val in values.items():
        rx = _var_re(name)
        n = len(rx.findall(text))
        if n != 1:
            sys.exit(f"atm_in: expected exactly one '{name}', found {n}")
        text = rx.sub(lambda m: m.group(1) + format_value(val) + m.group(3),
                      text)
    return text


# ------------------------------------------------------- nuopc.runconfig
def patch_runconfig(text, base_rundir, new_rundir, stop=None):
    """Point diro lines at the member dir; optionally override run length."""
    base = str(base_rundir).rstrip("/")
    n = text.count(base)
    text = text.replace(base, str(new_rundir))
    if stop:
        for key in ("stop_option", "stop_n"):
            rx = re.compile(rf"^(\s*{key}\s*=\s*)\S+", re.MULTILINE)
            if len(rx.findall(text)) != 1:
                sys.exit(f"nuopc.runconfig: expected one '{key}'")
            text = rx.sub(lambda m: m.group(1) + str(stop[key]), text)
    return text, n


# ------------------------------------------------------------ members.csv
def write_members(path, names, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["member"] + names)
        for i, row in enumerate(rows):
            w.writerow([member_name(i)] + [format_value(v) for v in row])


def read_members(path):
    with open(path, newline="") as f:
        r = csv.reader(f)
        header = next(r)
        names = header[1:]
        members = []
        for row in r:
            vals = {n: fortran_to_float(v) for n, v in zip(names, row[1:])}
            members.append((row[0], vals))
    return names, members
