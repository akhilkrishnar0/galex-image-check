import tempfile
import unittest
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from galex_image_check import CheckConfig,check_coverage,run_check,build_exclusions

class CoverageTests(unittest.TestCase):
    def setUp(self):
        self.w=WCS(naxis=2)
        self.w.wcs.crpix=[101,101]
        self.w.wcs.cdelt=[-1/60,1/60]
        self.w.wcs.crval=[150,20]
        self.w.wcs.ctype=['RA---TAN','DEC--TAN']
        self.data=np.ones((201,201))
        self.config=CheckConfig(radius_arcmin=10)
    def test_classifications(self):
        self.assertEqual(check_coverage(self.data,self.w,150,20,self.config)[0],'OK')
        ra,dec=self.w.pixel_to_world_values(1,100)
        self.assertEqual(check_coverage(self.data,self.w,ra,dec,self.config)[0],'EDGE')
        ra,dec=self.w.pixel_to_world_values(-30,100)
        self.assertEqual(check_coverage(self.data,self.w,ra,dec,self.config)[0],'OUTSIDE')
        self.data[100,100]=0
        self.assertEqual(check_coverage(self.data,self.w,150,20,self.config)[0],'INVALID')
        self.assertEqual(check_coverage(self.data,self.w,150,20,CheckConfig(10,3,.95,False))[0],'OK')
        self.data[100,100]=np.nan
        self.assertEqual(check_coverage(self.data,self.w,150,20,self.config)[0],'INVALID')
        self.assertEqual(check_coverage(self.data,self.w,np.nan,20,self.config)[0],'INVALID')
    def test_reports_and_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            folder=root/'inputs'/'TEST'/'FUV'/'intbgsub'
            folder.mkdir(parents=True)
            fits.writeto(folder/'good-fd-intbgsub.fits.gz',self.data,self.w.to_header())
            self.data[100,100]=0
            fits.writeto(folder/'bad-fd-intbgsub.fits',self.data,self.w.to_header())
            (folder/'broken-fd-intbgsub.fits').write_text('broken')
            cat=pd.DataFrame([('TEST',150,20)],columns=['galaxy','ra_deg','dec_deg'])
            with warnings.catch_warnings(record=True) as caught:
                table=run_check(cat,root/'inputs',root/'reports',self.config)
            self.assertTrue(caught)
            self.assertEqual(set(table.status),{'OK','INVALID','ERROR'})
            self.assertEqual(build_exclusions(table),{'TEST':['bad-fd-int.fits','broken-fd-int.fits']})
            for name in ['GALEX_GOOD_images.pdf','GALEX_INVALID_images.pdf','exclusions.json','settings.json','missing_observations.csv']:
                self.assertTrue((root/'reports'/name).is_file())
            with warnings.catch_warnings(record=True):
                empty=run_check(cat,root/'absent',root/'empty',make_pdfs=False)
            self.assertTrue(empty.empty)
            self.assertEqual(build_exclusions(empty),{})

if __name__=='__main__':
    unittest.main()
