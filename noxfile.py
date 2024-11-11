import nox
from laminci.nox import SYSTEM, build_docs, run, run_pre_commit, run_pytest

nox.options.default_venv_backend = "none"


@nox.session
def lint(session: nox.Session) -> None:
    run_pre_commit(session)


@nox.session
def build(session: nox.Session) -> None:
    run(session, f"uv pip install {SYSTEM} .[dev]")
    run_pytest(session, coverage=False)


# Currently not enabled
@nox.session
def docs(session: nox.Session) -> None:
    run(session, "lamin init --storage ./docsbuild --schema bionty,wetlab")
    build_docs(session, strict=True)
