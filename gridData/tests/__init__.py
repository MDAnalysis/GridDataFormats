# utilities for running the tests

import importlib


def module_not_found(module):
    try:
        importlib.import_module(module)
    except ImportError:
        return True
    else:
        return False
