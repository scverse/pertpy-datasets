import lamindb as ln
import pytest


@pytest.fixture(scope="session")
def lamindb_test_instance():
    ln.setup.init(storage="./test-pertpy-datasets", schema="bionty,wetlab,findrefs")
    yield
    ln.setup.delete("test-pertpy-datasets", force=True)
