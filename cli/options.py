"""Adapt Click's keyword arguments to typed command options."""

from functools import wraps
from typing import Callable, TypeVar

Options = TypeVar("Options")
Result = TypeVar("Result")


def command_options(options_type: type[Options]):
    """Place below Click option decorators to construct one options object."""

    def decorate(callback: Callable[[Options], Result]) -> Callable[..., Result]:
        @wraps(callback)
        def invoke(**kwargs) -> Result:
            return callback(options_type(**kwargs))

        return invoke

    return decorate
