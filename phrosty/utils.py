# IMPORTS Standard:
import os
import os.path as pa
import numpy as np
import pandas as pd
import requests
import warnings
from glob import glob

# IMPORTS Astro:
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
import astropy.wcs.wcs
from astropy.wcs.utils import skycoord_to_pixel
from astropy.table import Table
from astropy import units as u

# IMPORTS snpit
from snpit_utils.config import Config
from snpit_utils.logger import SNLogger

# The FITSFixedWarning is consequenceless and it will get thrown every single time you deal with a WCS.
warnings.simplefilter('ignore', category=FITSFixedWarning)


def _build_filepath(path, band, pointing, sca, filetype, rootdir=None):
    """Builds the filepath to an OpenUniverse simulation file.

    Parameters
    ----------
      path: Path
        If the path to the file is already known, this overrides the rest of the function and returns
        the input to kwarg 'path'.
      band: str
        Filter associated with target file.
      pointing: str
        Pointing number associated with target file.
      sca: str
        SCA number associated with target file.
      filetype: str
        The type of target file within the OpenUniverse simulations that you are looking for. Valid
        values are 'image' (*.fits.gz), 'truth' (*.fits.gz), and 'truthtxt' (*.txt).
      rootdir: Path, default None
        Root directory where OpenUniverse files are stored.

    Returns
    -------
      path: Path
        Path to target file.

    Raises
    ------
      ValueError
        if filetype is not 'image', 'truth', or 'truthtxt', a ValueError is raised.
      ValueError
        if (band is None) or (pointing is None) or (sca is None),
        you need to specify band, pointing, and sca if you do not provide a full filepath.
      ValueError
        if (path is None) and (band is None) and (pointing is None) and (sca is None),
        you need to provide either the full image path, or the band, pointing, and SCA.

    """

    rootdir = Config.get().value( 'ou24.tds_base' ) if rootdir is None else rootdir

    # First, what kind of file are we looking for?
    neededtypes = [ 'image', 'truth', 'truthtxt' ]
    if filetype not in neededtypes:
        raise ValueError(f"filetype must be in {neededtypes}.")
    elif filetype == 'image':
        subdir = 'images/simple_model'
        prefix = 'simple_model'
        extension = 'fits.gz'
    elif filetype == 'truth':
        subdir = 'images/truth'
        prefix = 'truth'
        extension = 'fits.gz'
    elif filetype == 'truthtxt':
        subdir = 'truth'
        prefix = 'index'
        extension = 'txt'

    # Did you already provide a path?
    if path is not None:
        return path
    elif (band is None) or (pointing is None) or (sca is None):
        raise ValueError('You need to specify band, pointing, and sca if you do not provide a full filepath.')
    elif (band is not None) and (pointing is not None) and (sca is not None):
        path = pa.join(rootdir, subdir, band, str(pointing),
                       f'Roman_TDS_{prefix}_{band}_{str(pointing)}_{str(sca)}.{extension}')
        return path

    elif (path is None) and (band is None) and (pointing is None) and (sca is None):
        raise ValueError('You need to provide either the full image path, or the band, pointing, and SCA.')


def ou2024_obseq_path( path=None ):
    """Retrieve the path to the OpenUniverse obseq file.

    Parameters
    ----------
      path: Path, default None
        Path to file. If not provided, use config file.

    Returns
    -------
    Path to OpenUniverse obseq file.

    """
    return ( os.path.join( Config.get().value('ou24.tds_base'), 'Roman_TDS_obseq_11_6_23.fits' )
             if path is None else path )


def get_roman_bands():
    """Get roman passbands.

    Returns
    -------
    List of bands included in the Roman-DESC TDS simulations.

    """
    return ['R062', 'Z087', 'Y106', 'J129', 'H158', 'F184', 'K213']


def read_truth_txt(path=None, band=None, pointing=None, sca=None):
    """Reads in the txt versions of the truth files as convenient astropy tables.

    Parameters
    ----------
      truthpath: str, default None
        Path to txt catalog version of truth file. If you do not
        provide this, you need to specify the arguments
        band, pointing, and sca.
      band: str, default None
        Roman filter. If you do not provide this, you need to provide truthpath.
      pointing: str, default None
        Pointing ID. If you do not provide this, you need to provide truthpath.
      sca: str, default None
        SCA ID. If you do not provide this, you need to provide truthpath.

    Returns
    -------
      truth: astropy.table.Table
        Astropy table with contents of the specified catalog txt file.

    """

    _truthpath = _build_filepath(path=path, band=band, pointing=pointing, sca=sca, filetype='truthtxt')
    truth_colnames = ['object_id', 'ra', 'dec', 'x', 'y', 'realized_flux', 'flux', 'mag', 'obj_type']
    truth_pd = pd.read_csv(_truthpath, comment='#', skipinitialspace=True, sep=' ', names=truth_colnames)
    truth = Table.from_pandas(truth_pd)

    return truth


def radec_isin(ra, dec, path=None, band=None, pointing=None, sca=None):
    """Check if a given RA, dec coordinate is in a given target image.

    Parameters
    ----------
    ra : float
        RA in degrees.
    dec : float
        Dec in degrees.
    path : Path, default None
        Path to image to check.
    band : str, default None
        Filter assocated with target image.
    pointing : str, default None
        Pointing associated with target image.
    sca : str, default None
        SCA associated with target image.

    Returns
    -------
    res: boolean
        True if provided RA, Dec is in the image. False if not.

    """

    _imgpath = _build_filepath(path, band, pointing, sca, 'image')
    with fits.open(_imgpath) as hdu:
        wcs = WCS(hdu[1].header)
    worldcoords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    try:
        x, y = skycoord_to_pixel(worldcoords, wcs)
    except astropy.wcs.wcs.NoConvergence:
        return False
    pxradec = np.array([x, y])
    if np.logical_or(any(pxradec < 0), any(pxradec > 4088)):
        res = False
    else:
        res = True

    return res


def get_corners(path=None, band=None, pointing=None, sca=None):
    """Retrieves the RA, dec of the corners of the specified SCA in degrees.

    :param band: Roman filter.
    :type band: str
    :param pointing: Pointing ID.
    :type pointing: str
    :param sca: SCA ID.
    :type sca: str
    :return: Tuple containing four numpy arrays, each with the RA and dec of the corner
            of the specified image in degrees.
    :rtype: tuple
    """
    _imgpath = _build_filepath(path, band, pointing, sca, 'image')
    with fits.open(_imgpath) as hdu:
        wcs = WCS(hdu[1].header)
    corners = [[0, 0], [0, 4088], [4088, 0], [4088, 4088]]
    wcs_corners = wcs.pixel_to_world_values(corners)

    return wcs_corners


# TODO clean this up for style
# TODO write something to clear this out
_parquet_cache = {}


def _read_parquet( file ):
    global _parquet_cache

    if file not in _parquet_cache:
        SNLogger.info( f"**** Reading parquet file {file}" )
        _parquet_cache[file] = pd.read_parquet( file )

        totm = 0
        for f, df in _parquet_cache.items():
            totm += df.memory_usage(index=True).sum()

        SNLogger.info( f"**** Done reading parquet file {file}; cache using {totm/1024/1024} MiB" )
    return _parquet_cache[ file ]


def get_transient_radec(oid):
    """Retrieve RA, dec of a transient based on its object ID."""

    oid = int(oid)
    snana_pq_path = os.path.join( Config.get().value('ou24.sn_truth_dir'), 'snana_*.parquet' )
    file_list = glob(snana_pq_path)
    for file in file_list:
        # Read the Parquet file
        df = _read_parquet(file)
        if len(df[df['id'] == oid]) != 0:
            ra = df[df['id'] == oid]['ra'].values[0]
            dec = df[df['id'] == oid]['dec'].values[0]
    return ra, dec


def get_transient_mjd(oid):
    """Retrieve start and end dates of a transient based on its object ID."""
    oid = int(oid)
    snana_pq_path = os.path.join( Config.get().value('ou24.sn_truth_dir'), 'snana_*.parquet' )
    file_list = glob(snana_pq_path)
    for file in file_list:
        # Read the Parquet file
        df = _read_parquet(file)
        if len(df[df['id'] == oid]) != 0:
            start = df[df['id'] == oid]['start_mjd'].values[0]
            end = df[df['id'] == oid]['end_mjd'].values[0]
    return start, end


def get_transient_zcmb(oid):
    """Retrieve z_CMB of a transient based on its object ID."""
    oid = int(oid)
    snana_pq_path = os.path.join( Config.get().value('ou24.sn_truth_dir'), 'snana_*.parquet' )
    file_list = glob(snana_pq_path)
    for file in file_list:
        # Read the Parquet file
        df = _read_parquet(file)
        if len(df[df['id'] == oid]) != 0:
            z = float(df[df['id'] == oid]['z_CMB'].values[0])

    return z


def get_transient_peakmjd(oid):
    """Retrieve z of a transient based on its object ID."""
    oid = int(oid)
    snana_pq_path = os.path.join( Config.get().value('ou24.sn_truth_dir'), 'snana_*.parquet' )
    file_list = glob(snana_pq_path)
    for file in file_list:
        # Read the Parquet file
        df = _read_parquet(file)
        if len(df[df['id'] == oid]) != 0:
            mjd = df[df['id'] == oid]['peak_mjd'].values[0]

    return mjd


def get_transient_info(oid):
    """Retrieve RA, Dec, MJD start, MJD end for specified object ID."""

    SNLogger.info( "*** calling get_transient_radec" )
    RA, DEC = get_transient_radec(oid)
    SNLogger.info( "*** calling get_transient_mjd" )
    start, end = get_transient_mjd(oid)
    SNLogger.info( "*** Done with get_transient_info" )

    return RA, DEC, start, end


def transient_in_or_out(oid, start, end, band):
    """Retrieve pointings that contain and do not contain the specified SN, per the truth files by MJD.

    transient_info_filepath is the output of make_object_table.

    Returns a tuple of astropy tables (images with the SN, images without the SN).
    """
    tab = Table.from_pandas( make_object_table( oid ) )

    tab.sort('pointing')
    tab = tab[tab['filter'] == band]

    in_all = get_mjd_info(start, end)
    in_rows = np.where(np.isin(tab['pointing'], in_all['pointing']))[0]
    in_tab = tab[in_rows]

    out_all = get_mjd_info(start, end, return_inverse=True)
    out_rows = np.where(np.isin(tab['pointing'], out_all['pointing']))[0]
    out_tab = tab[out_rows]

    return in_tab, out_tab


def get_mjd_limits(obseq_path=None):
    """Retrive the earliest and latest MJD in the simulations."""

    with fits.open(ou2024_obseq_path(obseq_path)) as obs:
        obseq = Table(obs[1].data)

    start = min(obseq['date'])
    end = max(obseq['date'])

    return start, end


def get_radec_limits(obseq_path=None):
    """Retrieve RA, dec limits (boresight coordinates)."""
    with fits.open(ou2024_obseq_path(obseq_path)) as obs:
        obseq = Table(obs[1].data)

    ra_min = min(obseq['ra'])
    ra_max = max(obseq['ra'])

    dec_min = min(obseq['dec'])
    dec_max = max(obseq['dec'])

    return {'ra': [ra_min, ra_max], 'dec': [dec_min, dec_max]}


def get_mjd(pointing, obseq_path=None):
    """Retrieve MJD of a given pointing.

    :param pointing: Pointing ID.
    :type pointing: int
    :param obseq_path: Path to obseq file Roman_TDS_obseq_11_6_23.fits.
    :type obseq_path: str, optional
    :return: MJD of specified pointing.
    :rtype: float
    """

    with fits.open(ou2024_obseq_path(obseq_path)) as obs:
        obseq = Table(obs[1].data)
    mjd = float(obseq['date'][int(pointing)])

    return mjd


def pointings_near_mjd(mjd, window=3, obseq_path=None):
    """Retrieve pointings near given MJD.

    :param mjd: Central MJD to search around.
    :type mjd: float
    :param window: Number of days around central MJD to include in search.
    :type window: float
    :param obseq_path: Path to obseq file Roman_TDS_obseq_11_6_23.fits.
    :type obseq_path: str, optional
    :return: Pointings within specified MJD range.
    :rtype: list
    """

    with fits.open(ou2024_obseq_path(obseq_path)) as obs:
        obseq = Table(obs[1].data)

    pointings = np.where(np.logical_and(obseq['date'] < mjd + window, obseq['date'] > mjd - window))[0]
    return pointings


def get_mjd_info(mjd_start=0, mjd_end=np.inf, return_inverse=False, obseq_path=None):
    """Get all pointings and corresponding filters between two MJDs.

    Returns an astropy table with columns 'filter' and 'pointing'.
    Does not return an 'sca' column because every sca belonging to a
    given pointing satisfies an MJD requirement.

    :param mjd_start: Start MJD, defaults to -np.inf
    :type mjd_start: float, optional
    :param mjd_end: End MJD, defaults to np.inf
    :type mjd_end: float, optional
    :param return_inverse: If true, returns all pointings outside the MJD range specified instead of inside.
    :type return_inverse: bool
    :param obseq_path: Path to obseq file Roman_TDS_obseq_11_6_23.fits.
    :type obseq_path: str, optional
    :return: Astropy table with pointing numbers and corresponding filters that satisfy the
            MJD requirements.
    :rtype: astropy.table.Table
    """
    with fits.open(ou2024_obseq_path(obseq_path)) as obs:
        obseq = Table(obs[1].data)

    if not return_inverse:
        mjd_idx = np.where((obseq['date'] > float(mjd_start)) & (obseq['date'] < float(mjd_end)))[0]
    elif return_inverse:
        mjd_idx = np.where(~((obseq['date'] > float(mjd_start)) & (obseq['date'] < float(mjd_end))))[0]

    mjd_tab = Table([obseq['filter'][mjd_idx], mjd_idx], names=('filter', 'pointing'))

    return mjd_tab


def get_exptime(band=None):

    exptime = {'F184': 901.175,
               'J129': 302.275,
               'H158': 302.275,
               'K213': 901.175,
               'R062': 161.025,
               'Y106': 302.275,
               'Z087': 101.7}

    if band in exptime.keys():
        return exptime[band]
    else:
        return exptime


def make_object_table(oid):

    ra, dec = get_transient_radec(oid)

    server_url = 'https://roman-desc-simdex.lbl.gov'
    req = requests.Session()
    result = req.post(f'{server_url}/findromanimages/containing=({ra}, {dec})')
    if result.status_code != 200:
        raise RuntimeError(f"Got status code {result.status_code}\n{result.text}")

    objs = pd.DataFrame(result.json())[['filter', 'pointing', 'sca']]
    return objs
