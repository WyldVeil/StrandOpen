"""Strand Open, Windows edition: lets the import pipeline (written for Linux)
run unchanged on Windows with the bundled tools. Python loads this file at
start (sitecustomize); on any other system it does nothing.

  - shell command lines (subprocess with shell=True) run in the bundled
    busybox shell, so pipes, quoting, 2>/dev/null, head and friends behave
    as on Linux; "bash" and "sh" scripts run there too
  - Unix tools the pipeline calls by name (sort, gzip, awk...) come from
    busybox; "hmmix" runs the bundled hmmix with this Python
  - text written to files and pipes uses "\\n" line ends, as on Linux (VCF
    lines, sample lists and region files must not gain a "\\r"); reading
    still accepts either
  - os.replace retries for a moment when the file is open elsewhere (the
    window reads the progress files the import rewrites)
"""
import os

if os.name == "nt":
    import builtins
    import io
    import shutil
    import subprocess
    import sys
    import time

    _BB = os.environ.get("STRANDOPEN_BUSYBOX") or os.path.join(os.path.dirname(sys.prefix), "tools", "busybox.exe")
    _APPLETS = {"sh", "bash", "awk", "sed", "sort", "cut", "head", "tail", "wc", "zcat", "gzip", "gunzip", "cat", "tr",
                "grep", "xargs", "tee", "seq", "uniq", "paste", "join", "comm", "tac", "mv", "rm", "mkdir", "ls", "touch",
                "nproc", "dirname", "basename", "sleep", "printf", "echo"}
    _HMMIX = "import sys; from main import main; sys.argv[0] = 'hmmix'; sys.exit(main())"

    # ---- commands
    _popen_init = subprocess.Popen.__init__

    def _init(self, args, *a, **kw):
        if kw.get("shell") and isinstance(args, str):
            args, kw["shell"] = [_BB, "sh", "-c", args], False
            kw.pop("executable", None)
        elif isinstance(args, (list, tuple)) and args:
            first = os.path.basename(str(args[0]))
            if first == "hmmix":
                args = [sys.executable, "-c", _HMMIX] + list(args[1:])
            elif first in ("bash", "sh"):
                args = [_BB, "sh"] + list(args[1:])
            elif first in _APPLETS and shutil.which(first) is None:
                args = [_BB, first] + list(args[1:])
            if kw.get("executable") in ("/bin/bash", "/bin/sh"):
                kw.pop("executable")
        _popen_init(self, args, *a, **kw)

    subprocess.Popen.__init__ = _init

    # ---- "\n" line ends when writing text
    _open = builtins.open

    def _lf_open(file, mode="r", buffering=-1, encoding=None, errors=None, newline=None, *a, **kw):
        if newline is None and "b" not in mode and "+" not in mode and any(c in mode for c in "wax"):
            newline = "\n"
        return _open(file, mode, buffering, encoding, errors, newline, *a, **kw)

    builtins.open = io.open = _lf_open

    class _LFTextIOWrapper(io.TextIOWrapper):
        """write-only text streams (pipes into a program, gzip.open "wt") write "\\n" ends"""
        def __init__(self, buffer, encoding=None, errors=None, newline=None, *a, **kw):
            try:
                write_only = buffer.writable() and not buffer.readable()
            except (AttributeError, ValueError, OSError):
                write_only = False
            if newline is None and write_only:
                newline = "\n"
            super().__init__(buffer, encoding, errors, newline, *a, **kw)

    io.TextIOWrapper = _LFTextIOWrapper

    # ---- replacing a file another process has open
    _replace = os.replace

    def _retry_replace(src, dst, *a, **kw):
        for i in range(40):
            try:
                return _replace(src, dst, *a, **kw)
            except PermissionError:
                if i == 39:
                    raise
                time.sleep(0.05 * (i + 1) if i < 10 else 0.5)

    os.replace = _retry_replace
