"""Stub package: force i-PI to skip PyFFTW and use NumPy FFT.

i-PI only catches ImportError around ``import pyfftw``; a broken PyFFTW install
can raise ValueError during allocation instead.  Raising ImportError here makes
``nmtransform.nm_fft`` fall back to NumPy (slower but correct).

The orchestrator prepends this directory to PYTHONPATH for the ``i-pi``
subprocess only when a quick probe detects a broken PyFFTW.
"""

raise ImportError(
    "pyfftw disabled by uma_ice_rpmd orchestrator probe (broken n_byte_align_empty); "
    "i-PI will use NumPy FFT."
)
