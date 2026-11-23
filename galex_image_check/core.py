"""Local GALEX field-coverage screening; never deletes observations."""
from pathlib import Path
from dataclasses import dataclass,asdict
import json
import warnings
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

@dataclass(frozen=True)
class CheckConfig:
    radius_arcmin: float=3
    plot_arcmin: float=3.0
    min_coverage: float=100
    zero_is_invalid: bool=True
    def __post_init__(self):
        if not np.isfinite([self.radius_arcmin,self.plot_arcmin,self.min_coverage]).all():
            raise ValueError('Settings must be finite')
        if self.radius_arcmin<=0 or self.plot_arcmin<=0 or not 0<=self.min_coverage<=1:
            raise ValueError('Require positive radii and coverage between 0 and 1')

def load_fits(path):
    with fits.open(path,memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is None or hdu.data.ndim!=2:
                continue
            wcs=WCS(hdu.header).celestial
            if wcs.has_celestial:
                return np.array(hdu.data,dtype=float,copy=True),wcs
    raise ValueError('No valid 2D celestial FITS image')

def get_scale(wcs):
    scale=float(np.mean(np.abs(proj_plane_pixel_scales(wcs)))*60)
    if not np.isfinite(scale) or scale<=0:
        raise ValueError('Invalid WCS pixel scale')
    return scale

def check_coverage(data,wcs,ra,dec,config=None):
    config=config or CheckConfig()
    ny,nx=data.shape
    x,y=map(float,wcs.world_to_pixel_values(ra,dec))
    if not np.isfinite([x,y]).all():
        return 'INVALID',x,y,0.0
    rp=config.radius_arcmin/get_scale(wcs)
    if not (0<=x<nx and 0<=y<ny):
        return 'OUTSIDE',x,y,0.0
    x1,x2=max(0,int(np.floor(x-rp-1))),min(nx,int(np.ceil(x+rp+2)))
    y1,y2=max(0,int(np.floor(y-rp-1))),min(ny,int(np.ceil(y+rp+2)))
    yy,xx=np.mgrid[y1:y2,x1:x2]
    region=(xx-x)**2+(yy-y)**2<=rp**2
    values=data[y1:y2,x1:x2]
    valid=np.isfinite(values)
    if config.zero_is_invalid:
        valid &= values!=0
    # Preserve the supplied notebook's continuous-area denominator.
    fraction=min(1.0,np.count_nonzero(valid&region)/max(1,int(np.ceil(np.pi*rp**2))))
    cx,cy=int(np.rint(x)),int(np.rint(y))
    if not (0<=cx<nx and 0<=cy<ny):
        return 'OUTSIDE',x,y,fraction
    if not np.isfinite(data[cy,cx]) or (config.zero_is_invalid and data[cy,cx]==0):
        return 'INVALID',x,y,fraction
    return ('EDGE' if fraction<config.min_coverage else 'OK'),x,y,fraction

def build_exclusions(table):
    """Convert flagged background-subtracted names to pipeline intensity basenames."""
    result={str(g):[] for g in table.galaxy.unique()}
    for row in table.loc[table.status.ne('OK')].itertuples():
        name=Path(row.file).name.removesuffix('.gz')
        if not name.endswith('-intbgsub.fits'):
            raise ValueError(f'Unexpected observation filename: {name}')
        result[str(row.galaxy)].append(name.replace('-intbgsub.fits','-int.fits'))
    return {g:sorted(set(names)) for g,names in result.items()}

def run_check(catalog,base_dir='downloaded_galex',output_dir='galex_image_check_results',config=None,make_pdfs=True):
    """Scan FUV/NUV intbgsub folders; return a table and write CSV/JSON/PDF reports."""
    config=config or CheckConfig()
    catalog=pd.DataFrame(catalog).copy()
    if catalog.empty or not {'galaxy','ra_deg','dec_deg'}<=set(catalog):
        raise ValueError('Supply a nonempty catalog with galaxy, ra_deg, dec_deg')
    coords=catalog[['ra_deg','dec_deg']].to_numpy(dtype=float)
    if not np.isfinite(coords).all() or not ((coords[:,0]>=0)&(coords[:,0]<360)&(abs(coords[:,1])<=90)).all():
        raise ValueError('Invalid RA/Dec in degrees')
    if catalog.galaxy.isna().any():
        raise ValueError('Missing galaxy name')
    catalog['galaxy']=catalog.galaxy.astype(str).str.replace(' ','',regex=False)
    if catalog.galaxy.duplicated().any() or any(not n or n in ('.','..') or '/' in n or '\\' in n for n in catalog.galaxy):
        raise ValueError('Galaxy names must be unique and file-safe')
    rows,missing=[],[]
    for gal in catalog.itertuples():
        for band in ('FUV','NUV'):
            folder=Path(base_dir)/gal.galaxy/band/'intbgsub'
            files=sorted(set(folder.glob('*-intbgsub.fits'))|set(folder.glob('*-intbgsub.fits.gz')))
            if not files:
                missing.append(dict(galaxy=gal.galaxy,band=band,folder=str(folder)))
                warnings.warn(f'Missing observations: {folder}',RuntimeWarning)
            for file in files:
                error=''
                try:
                    data,wcs=load_fits(file)
                    status,x,y,fraction=check_coverage(data,wcs,gal.ra_deg,gal.dec_deg,config)
                except Exception as exc:
                    status,x,y,fraction='ERROR',np.nan,np.nan,0.0
                    error=f'{type(exc).__name__}: {exc}'
                rows.append(dict(galaxy=gal.galaxy,band=band,file=str(file),ra_deg=gal.ra_deg,dec_deg=gal.dec_deg,x_pixel=x,y_pixel=y,coverage_fraction=fraction,status=status,error=error))
    table=pd.DataFrame(rows,columns=['galaxy','band','file','ra_deg','dec_deg','x_pixel','y_pixel','coverage_fraction','status','error'])
    out=Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    table.to_csv(out/'GALEX_field_edge_check.csv',index=False)
    pd.DataFrame(missing,columns=['galaxy','band','folder']).to_csv(out/'missing_observations.csv',index=False)
    (out/'settings.json').write_text(json.dumps(asdict(config),indent=2)+'\n')
    (out/'exclusions.json').write_text(json.dumps(build_exclusions(table),indent=2)+'\n')
    if make_pdfs:
        from .plotting import create_pdf
        create_pdf(table,['OK'],out/'GALEX_GOOD_images.pdf',config)
        create_pdf(table,['EDGE','OUTSIDE','INVALID','ERROR'],out/'GALEX_INVALID_images.pdf',config)
    return table
