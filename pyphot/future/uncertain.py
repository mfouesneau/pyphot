"""Implements what could be a filter with uncertainties

This is not a complete implementation yet. It is basically a Monte-Carlo simulation of a filter with uncertainties.
"""

from typing import Literal, Union, Optional, Any, Sequence, Tuple

import numpy as np
import numpy.typing as npt
from scipy.integrate import trapezoid

from . import config
from .unit_adapters import QuantityType


__all__ = ["Filter"]
from .phot import Filter, _drop_units, _split_value_unit


class UncertainFilter(Filter):
    """Implements what could be a filter with uncertainties

    Attributes
    ----------
    wavelength: ndarray
        wavelength sequence defining the filter transmission curve

    mean: Filter
        mean passband transmission

    samples: sequence(Filter)
        samples from the uncertain passband transmission model

    name: string
        name of the passband

    dtype: str
        detector type, either "photon" or "energy" counter

    unit: str
        wavelength units

    vega: str
        Vega flavor to use for calculations, default is 'default'
        (see :class:`pyphot.vega.Vega` for details)
    """

    dtype: Literal["photon", "energy"]
    name: str
    norm: float
    transmit: npt.NDArray[np.floating]
    wavelength_unit: str
    samples_: Sequence[Filter]
    model_: Optional[Any] = None

    _cl: float
    _lT: float
    _lpivot: float
    _vega_flavor: str
    _wavelength: npt.NDArray[np.floating]

    def __init__(
        self,
        wavelength: Union[npt.NDArray[np.floating], QuantityType],
        mean_transmit: npt.NDArray[np.floating],
        sample_transmit: Sequence[npt.NDArray[np.floating]],
        *,
        name: str = "",
        dtype: Literal["photon", "energy"] = "photon",
        unit: Optional[str] = None,
        vega: Optional[str] = None,
    ):
        """Constructor"""
        super().__init__(wavelength, mean_transmit, name=name, dtype=dtype, unit=unit)
        self.samples_ = [
            Filter(
                wavelength,
                transmit_k,
                name=name + "_{0:d}".format(num),
                dtype=dtype,
                unit=unit,
                vega=vega,
            )
            for (num, transmit_k) in enumerate(sample_transmit)
        ]
        self.model_ = None

    @classmethod
    def from_gp_model(
        cls,
        model,
        xprime=None,
        n_samples: int = 10,
        **kwargs,
    ) -> "UncertainFilter":
        """Generate a filter object from a sklearn GP model

        Parameters
        ----------
        model: sklearn.gaussian_process.GaussianProcessRegressor
            model of the passband
        xprime: ndarray
            wavelength to express the model in addition to the training points
        n_samples: int
            number of samples to generate from the model.
        kwargs: dict
            UncertainFilter keywords
        """
        if xprime is None:
            xpred = model.X_train_
        else:
            xpred = np.unique(np.hstack([_drop_units(xprime), model.X_train_.ravel()]))
            xpred = xpred.reshape(1, -1).T

        unit_ = kwargs.pop("unit", None)
        if unit_ is None:
            unit_ = str(getattr(xprime, "units", None))

        mean_transmit, _ = model.predict(xpred, return_std=True)
        samples = model.sample_y(xpred, n_samples=n_samples)

        unc_filter = cls(xpred.ravel(), mean_transmit, samples.T, unit=unit_, **kwargs)

        unc_filter.model_ = model
        return unc_filter

    def info(self, show_zeropoints=True):
        """display information about the current filter"""
        string = self._get_info(show_zeropoints)
        string = string.replace(
            "Filter object information", "Filter object mean information only"
        )
        print(string)

    def _get_samples_attribute(
        self, attr: str, *args, **kwargs
    ) -> Union[npt.NDArray[Any], QuantityType]:
        """Returns the attribute from all samples"""
        try:  # callable attribute
            vals = [getattr(fk, attr)(*args, **kwargs) for fk in self.samples_]
        except TypeError:
            vals = [getattr(fk, attr) for fk in self.samples_]
        _, unit_ = _split_value_unit(vals[0])
        if unit_ is not None:
            unit_ = config.units.U(unit_)
            return np.array([_drop_units(v) for v in vals]) * unit_
        else:
            return np.array(vals)

    def _get_mean_attribute(
        self, attr: str, *args, **kwargs
    ) -> Union[npt.NDArray[Any], QuantityType]:
        """Returns the attribute from the mean passband"""
        attr = getattr(self, attr)
        try:
            return attr(*args, **kwargs)  # type: ignore / try
        except TypeError:
            return attr

    def _get_mean_and_samples_attribute(
        self, attr: str, *args, **kwargs
    ) -> Tuple[
        Union[npt.NDArray[Any], QuantityType],
        Union[npt.NDArray[Any], QuantityType],
    ]:
        """Compute / extract mean and smapled filter attributes

        Parameters
        ----------
        attr: str
            attribute to get (can be a callable attribute)
        args: sequence
            any argument of attr
        kwargs: dict
            any keywords for attr

        Returns
        -------
        mean_: object
            value from the mean passband
        samples_: sequence(object)
            values from each sampled passband
        """
        return (
            self._get_mean_attribute(attr, *args, **kwargs),
            self._get_samples_attribute(attr, *args, **kwargs),
        )

    @property
    def transmit_mean_and_samples(self):
        """Transmission curves"""
        return self._get_mean_and_samples_attribute("transmit")

    def reinterp(
        self, lamb: Union[npt.NDArray[np.floating], QuantityType]
    ) -> "UncertainFilter":
        """reinterpolate filter onto a different wavelength definition

        Parameters
        ----------
        lamb : QuantityType | npt.NDArray[np.floating]
            Wavelength to reinterpolate onto
            If no unit is provided, the filter's wavelength unit is assumed.

        Returns
        -------
        UncertainFilter
            Filter reinterpolated onto the new wavelength definition
        """
        mean = Filter.reinterp(self, lamb)
        samples = [fk.reinterp(lamb) for fk in self.samples_]

        wavelengths = mean.wavelength
        transmissions = mean.transmit
        sample_transmissions = [sk.transmit for sk in samples]

        return self.__class__(
            wavelengths,
            transmissions,
            sample_transmissions,
            name=self.name,
            dtype=self.dtype,
            unit=mean.wavelength_unit,
            vega=self._vega_flavor,
        )
