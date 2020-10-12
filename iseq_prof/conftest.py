def pytest_addoption(parser):
    try:
        parser.addoption(
            "--runslow",
            action="store_true",
            default=False,
            help="run slow tests",
        )
    except ValueError as e:
        if "already added" not in e.args[0]:
            raise e


def pytest_collection_modifyitems(config, items):
    import pytest

    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return

    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
