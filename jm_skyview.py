#!/usr/bin/python
"""
Script to retrieve cutout images from NASA Skyview service from fields identified by source name or coordinates.

Some of the functions are adapted from https://github.com/cosmicpudding/skyviewbot/tree/master/skyviewbot

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
"""

import glob
import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
import astropy.units as u
import aplpy
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy

def coords_from_name(field_name):
    """Get ra, dec coordinates from a field name using astropy

    Args:
        field_name (str): Field name, e.g. 'M101'

    Returns:
        (float, float): ra, dec in degrees

    Example:
        >>> coords_from_name('M101')
        (210.80242917, 54.34875)
    """
    coord = SkyCoord.from_name(field_name)

    return coord.ra.to(u.deg).value, coord.dec.to(u.deg).value


def plot_fits(fits_name, plot_title=None, cmap_name='viridis', colorbar=True, contour=True):
    """Make a PNG plot out of a FITS file
    
    Args:
        fits_name (str): path of fits file
        plot_title (str): plot title, default is name of the fits file
        cmap_name (str): name of colormap, default is viridis
        colorbar (bool): include colorbar, default is True
        contour (bool): include contour, default is True
    """
    f = aplpy.FITSFigure(fits_name, figsize=(10, 8))
    if plot_title == None:
        plot_title = fits_name.replace('.fits', '')
    plt.title(plot_title)
    f.show_colorscale(cmap=cmap_name, stretch='linear')
    f.ticks.set_color('k')
    if colorbar:
        f.add_colorbar()
    if 'BMAJ' in fits.open(fits_name)[0].header:
        f.add_beam()
        print(f'Adding beam for {fits_name}')
    if contour:
        f.show_contour()
    output_name = fits_name.replace('.fits', '.png')
    plt.savefig(output_name, dpi=200, bbox_inches='tight')

def call_skyview(survey, pos, fov, coord, fitsname, proj='Car', pix=500):
    """Call Skyview to download data from a survey based on input parameters

    Args:
        survey (str): name of survey, from https://skyview.gsfc.nasa.gov/current/cgi/survey.pl
        pos (float,float): position coordinates as a tuple
        fov (float): FOV in degrees
        coord (str): coordinate system (e.g. Galactic, J2000, B1950)
        fitsname (str): name of output fits file
        proj (str): projection of image. (e.g. Car, Sin)
        pix (int): pixel dimensions of image (e.g. 500)

    Examples:
        >>> call_skyview('DSS', (255.291,-29.911), 5., 'J2000', '/tmp/bla.fits')
        >>> call_skyview('NVSS', (144.99497,-01.22029), 0.5, 'Gal', '/tmp/bla.fits')
    """

    x, y = pos

    images = SkyView.get_images(SkyCoord(ra=x*u.deg, dec=y*u.deg), survey,
                                coordinates=coord,
                                projection=proj, pixels=pix,
                                height=fov*u.deg, width=fov*u.deg)

    try:
        images[0][0].writeto(fitsname, overwrite=True)
    except astropy.io.fits.verify.VerifyError:
        print('Data not available')
        pass
    return fitsname 

def main(*function_args):
    """Command-line interface to skyviewbot

    Args:
        List[str]: arguments (typically sys.argv[1:])

    Returns:
        bool: True if succeeded, False if not
    """

    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('field', help='Field, e.g. "M101" or "255.2,1" (if it contains a comma, '
                                      'it\'s interpreted as coordinates, otherwise fed to CDS)')
    parser.add_argument('-s', '--survey',
                        default='DSS',
                        type=str,
                        help="Survey name, e.g. 'DSS' (default: %(default)s)")
    parser.add_argument('-r', '--radius',
                        default=1.,
                        type=float,
                        help="Radius (default: %(default)s)")
    parser.add_argument('-c', '--colormap',
                        default="viridis",
                        type=str,
                        help="Colormap (default: %(default)s)")
    parser.add_argument('-p', '--path',
                        default="./",
                        type=str,
                        help="path to save images (default: %(default)s)")
    args = parser.parse_args(*function_args)


    field = args.field
    print(f'\nFIELD: {field}')
    print(f'SURVEY: {args.survey}')
    if ',' in field:
        ra_str, dec_str = field.split(',')
        ra, dec = float(ra_str), float(dec_str)
        fieldname = f'ra{ra}dec{dec}'
    else:
        ra, dec = coords_from_name(field)
        fieldname = field

    survey = args.survey.replace(' ', '_')
    fits_name = f'{fieldname}_{survey}_{args.radius}d.fits'
    fits_path = os.path.join(args.path, fits_name)
    if not os.path.isfile(fits_path):
        print('Downloading file')
        call_skyview(args.survey, (ra, dec), args.radius, 'J2000', fits_path)
    else:
        print(f'Already downlaoded images: {fits_path}')
    if os.path.isfile(fits_path) and not os.path.isfile(fits_path.replace('.fits', '.png')):
        plot_fits(fits_path, plot_title=None, cmap_name='viridis', colorbar=True)


if __name__ == '__main__':
    main(sys.argv[1:])

