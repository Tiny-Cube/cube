"""Electronic density of states with explicit state and spin normalization."""
__all__ = ['broaden', 'kmesh', 'bloch_spectra', 'bloch_eigensystems',
           'participation_groups', 'disorder_statistics']


def __getattr__(name):
    if name in __all__:
        from . import core
        return getattr(core, name)
    raise AttributeError(name)
