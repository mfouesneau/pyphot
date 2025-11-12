"""
UnitsAdapter base class and decorator
=====================================

Provides the common interface for units adapters.
"""

from dataclasses import dataclass
from typing import Any, Optional, Callable, Sequence, cast, Union
from functools import wraps
import warnings
import inspect
from textwrap import dedent


__all__ = ["UnitsAdapter", "UnitTyping"]


def _warning_on_one_line(
    message: str, category: Any, filename: str, lineno: int, file=None, line=None
) -> str:
    """Prints a complete warning that includes exactly the code line triggering it from the stack trace."""
    return " {:s}:{:d} {:s}:{:s}".format(
        filename, lineno, category.__name__, str(message)
    )


def raise_warning(
    msg: str,
    formatwarning: Optional[Callable] = _warning_on_one_line,
    **kwargs,
):
    """Raise a warning with a custom format local format.

    Parameters
    ----------
    msg: str
        The warning message.
    formatwarning: Callable, optional
        The function to format the warning message. Defaults to _warning_on_one_line.
    **kwargs:
        Additional keyword arguments to pass to the warning function.
    """
    if formatwarning is None:
        warnings.warn(msg, **kwargs)
    else:
        prev_formatwarning = warnings.formatwarning
        warnings.formatwarning = _warning_on_one_line
        warnings.warn(msg, **kwargs)
        warnings.formatwarning = prev_formatwarning


@dataclass
class UnitTyping:
    """Collects exceptions related to the unit handling as standard exceptions"""

    DimensionalityError: Any
    UndefinedUnitError: Any
    Quantity: Any


class UnitsAdapter:
    """Unifies unit handling"""

    typing: UnitTyping

    @staticmethod
    def U(*args, **kwargs) -> Any:
        """Returns the resulting object from a unit query to the unit registry.

        Note:
            Some registries make a distinction between units and quantities (i.e. value with unit information).

        seealso: UnitsAdapter.Q
        """
        raise NotImplementedError

    @staticmethod
    def Q(*args, **kwargs) -> Any:
        """Returns an explicit quantity from a unit query (not always the default)"""
        raise NotImplementedError

    @staticmethod
    def has_unit(val: Any) -> bool:
        """Check if val is a value with unit information"""
        return hasattr(val, "units") or hasattr(val, "unit")

    @staticmethod
    def val_in_unit(
        varname: str,
        value: Any,
        defaultunit: Optional[str] = None,
        warn: bool = True,
    ) -> Any:
        raise NotImplementedError

    @classmethod
    def decorate(
        cls,
        *args: str,
        output: Optional[str] = None,
        warn: bool = False,
        **kwargs,
    ) -> Callable:
        """Decorator to enforce default units for function arguments.

        Parameters
        ----------
        *args : str
            Default units for function arguments.
        output : str, optional
            Default unit for the function output.
        warn : bool, optional
            Whether to warn if a variable does not have explicit units.
        _Adapter : UnitsAdapter, optional
            Units adapter to impose, defaulting to the global adapter.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        Callable
            Decorated function.
        """
        return enforce_default_units(
            *args,
            output=output,
            warn=warn,
            Adapter=cast(UnitsAdapter, cls),
            **kwargs,
        )


class enforce_default_units:
    def __init__(
        self,
        *args: Union[str, None],
        output: Optional[str] = None,
        warn: bool = False,
        Adapter: Optional[UnitsAdapter] = None,
        **kwargs,
    ):
        """Decorator to enforce default units for function arguments.

        .. note::
            This decorator updates the docstring of the decorated function to reflect the enforced units.

        Parameters
        ----------
        *args : str
            Default units for function arguments.
        output : str, optional
            Default unit for the function output.
        warn : bool, optional
            Whether to warn if a variable does not have explicit units.
        Adapter : UnitsAdapter, optional
            Units adapter to impose, defaulting to the global adapter.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        Callable
            Decorated function.
        """
        self._adapter: Optional[UnitsAdapter] = Adapter
        self.warn: bool = warn
        self.args: Sequence[Union[str, None]] = args
        self.kwargs: dict = kwargs
        self.output: Optional[str] = output
        self.fn_signature: Optional[inspect.Signature] = None
        self.fn_argnames: Optional[Sequence[str]] = None
        self.func: Optional[Callable] = None

    @property
    def adapter(self) -> UnitsAdapter:
        # I am not a fan of this hack, but it's the only way to avoid circular imports
        from .. import config

        return self._adapter or config.units

    def parse_function_arguments(self, fn_args, fn_kwargs):
        if self.func is None:
            raise ValueError("Function not set")
        if (self.fn_argnames is None) or (self.fn_signature is None):
            self.fn_signature = inspect.signature(self.func)
            self.fn_argnames = list(self.fn_signature.parameters.keys())

        adapter = self.adapter
        fn_argnames = self.fn_argnames

        n_args = min(len(fn_args), len(self.args))
        new_args = []
        # assume args correspond to the first arguments of the function
        for arg_name, arg_value, unit in zip(fn_argnames[:n_args], fn_args, self.args):
            if arg_value is None:
                continue
            if isinstance(arg_value, str) and arg_value.strip() == "":
                continue
            new_args.append(
                adapter.val_in_unit(
                    arg_name,
                    arg_value,
                    unit,
                    warn=self.warn,
                )
            )
        # add other arguments
        new_args.extend(fn_args[n_args:])

        # deal with keywords
        new_kwargs = {
            par.name: par.default
            for par in self.fn_signature.parameters.values()
            if par.default is not inspect._empty
        }
        new_kwargs.update(fn_kwargs)
        for key, unit in self.kwargs.items():
            new_kwargs[key] = adapter.val_in_unit(
                key, new_kwargs[key], unit, warn=self.warn
            )

        return new_args, new_kwargs

    def generate_docstring(self) -> str:
        """Generate a docstring for the decorated function."""

        if self.func is None:
            raise ValueError("Function not set")
        if (self.fn_argnames is None) or (self.fn_signature is None):
            self.fn_signature = inspect.signature(self.func)
            self.fn_argnames = list(self.fn_signature.parameters.keys())

        docstring = []
        if self.args or self.kwargs or self.output:
            docstring.append(""".. important::""")
            docstring.append(
                f"   Decorated function :func:`{self.func.__name__}` with default units.\n"
            )
        if self.args:
            docstring.append("   Parameters Units:\n")
            for name, unit in zip(self.fn_argnames, self.args):
                if unit in (None, ""):
                    continue
                docstring.append(f"   - {name} : {unit}")
            docstring.append("")
        if self.kwargs:
            docstring.append("   Keywords Units\n")
            for name, unit in self.kwargs.items():
                docstring.append(f"   - {name} : {unit}")
            docstring.append("")
        if self.output:
            docstring.append("   Return Units\n")
            docstring.append(f"   - output: {self.output}")
            docstring.append("")
        __deco_doc__ = "\n".join(docstring)

        if self.func.__doc__:
            fn_doc = dedent(self.func.__doc__.strip()).splitlines()
            __doc__ = "\n".join(
                [
                    fn_doc[0].strip(),
                    "",
                    dedent(__deco_doc__.strip()),
                    dedent("\n".join(fn_doc[1:])),
                ]
            )
        else:
            __doc__ = __deco_doc__
        return __doc__

    def __call__(self, func: Callable) -> Callable:
        """
        Decorate a function to enforce default units.

        Parameters
        ----------
        func : Callable
            Function to decorate.

        Returns
        -------
        Callable
            Decorated function.
        """
        self.func = func

        @wraps(func)
        def wrapper(*fn_args, **fn_kwargs):
            # Enforce default units for function arguments
            new_args, new_kwargs = self.parse_function_arguments(fn_args, fn_kwargs)

            # Call the function with the adapted arguments
            result = func(*new_args, **new_kwargs)

            # Force the output to the desired unit
            return self.adapter.val_in_unit(
                f"{func.__name__}(...)", result, self.output
            )

        # update the docstring accordingly
        wrapper.__doc__ = self.generate_docstring()

        return wrapper
