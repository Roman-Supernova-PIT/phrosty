__all__ = [ 'band_aliases', 'map_band', 'managed_pool' ]

import contextlib
import gc
from multiprocessing import Pool

import nvtx

from snappl.logger import SNLogger

_band_dict = {
                "F062": "F062",
                "F087": "F087",
                "F106": "F106",
                "F129": "F129",
                "F158": "F158",
                "F184": "F184",
                "F213": "F213",
                "R062": "F062",
                "Z087": "F087",
                "Y106": "F106",
                "J129": "F129",
                "H158": "F158",
                "K213": "F213",
            }


def map_band(band_input):
    """Maps the possible band inputs in two different
    formats for flexibility.
    """

    try:
        return _band_dict[band_input]
    except KeyError:
        raise ValueError(f'Unknown band {band_input}; must be one of {list(_band_dict)}.')


def band_aliases( band_input ):
    """Return every accepted spelling of a band, with band_input first.

    e.g. band_aliases('F106') -> ['F106', 'Y106']
         band_aliases('Y106') -> ['Y106', 'F106']

    Written by Claude Opus 5.5.

    """
    canonical = map_band( band_input )
    aliases = [ b for b, f in _band_dict.items() if f == canonical ]
    # Try the spelling the user actually wrote first
    aliases.sort( key=lambda b: b != band_input )
    return aliases


@contextlib.contextmanager
def managed_pool( nprocs, wait_label=None, wait_color=0xff8888 ):
    """Context manager for a multiprocessing Pool that cleans up deterministically.

    This function was written by Claude Opus 5.5.

    Parameters
    ----------
      nprocs : int
        Number of worker processes. If nprocs <= 1, no pool is created and
        None is yielded, so callers can use the same `with` statement for
        both serial and parallel runs.

      wait_label : str or None, default None
        If given, waiting for the workers to finish (close + join) is wrapped
        in an nvtx range with this name, and logged. Use it to profile how
        long the pool's outstanding work takes after you stop submitting.

      wait_color : int, default 0xff8888
        Color of the nvtx range.

    Yields
    ------
      multiprocessing.pool.Pool or None

    Notes
    -----
    Pool objects sit in reference cycles, so their pipes can stay open until
    the garbage collector happens to run. Over many pools this exhausts the
    file descriptor limit (256 by default on macOS). Collecting explicitly
    on exit prevents that.

    You do not need to call close() or join() inside the `with` block;
    leaving the block does it for you.
    """

    # Serial mode: nothing to create or clean up
    if nprocs <= 1:
        yield None
        return

    pool = Pool( nprocs )
    try:
        yield pool

    except BaseException:
        # The body failed; don't wait on queued work, just kill the workers
        pool.terminate()
        pool.join()
        raise

    else:
        # Normal exit: stop accepting work, then wait for outstanding tasks.
        # This wait is the part worth profiling.
        wait_ctx = ( nvtx.annotate( wait_label, color=wait_color )
                     if wait_label is not None else contextlib.nullcontext() )
        with wait_ctx:
            if wait_label is not None:
                SNLogger.info( f"Waiting for {wait_label} processes to finish..." )
            pool.close()
            pool.join()
            if wait_label is not None:
                SNLogger.info( f"...{wait_label} processes done." )

        # Tear down the pool's handler threads
        pool.terminate()

    finally:
        # Break reference cycles so the pool's pipes are actually closed now
        del pool
        gc.collect()
