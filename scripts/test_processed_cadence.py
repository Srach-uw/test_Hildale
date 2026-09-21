import numpy as np
import pytest
from astropy.io import fits
from verify_processed_cadence import inspect, verify_products


def make(path,seconds,errors=None):
    fits.HDUList([fits.PrimaryHDU(),
        fits.ImageHDU(np.arange(20)*seconds/86400,name='TIME'),
        fits.ImageHDU(np.ones(20),name='FLUX'),
        fits.ImageHDU(np.ones(20)*.001 if errors is None else errors,name='ERROR')]).writeto(path)


def test_short_cadence_measurement(tmp_path):
    path=tmp_path/'sc.fits'; make(path,58.85)
    assert inspect(path,'short')['usable_rows']==20


def test_renaming_long_as_short_is_rejected(tmp_path):
    path=tmp_path/'sc.fits'; make(path,1765.5)
    with pytest.raises(ValueError,match='cadence mismatch'):
        inspect(path,'short')


def test_invalid_errors_rejected(tmp_path):
    path=tmp_path/'sc.fits'; make(path,58.85,np.zeros(20))
    with pytest.raises(ValueError,match='Insufficient'):
        inspect(path,'short')


def test_detrending_alone_does_not_validate_sampler_product(tmp_path):
    make(tmp_path/'K00001_sc_detrended.fits',58.85)
    with pytest.raises(ValueError,match='filtered'):
        verify_products(tmp_path,'K00001','both',['detrended','filtered'])


def test_both_stages_verified(tmp_path):
    for stage in ['detrended','filtered']:
        make(tmp_path/f'K00001_sc_{stage}.fits',58.85)
    result=verify_products(tmp_path,'K00001','both',['detrended','filtered'])
    assert result['filtered']['short']['usable_rows']==20


def test_stale_sc_rejected_in_lc_control(tmp_path):
    make(tmp_path/'K00001_lc_filtered.fits',1765.5)
    make(tmp_path/'K00001_sc_filtered.fits',58.85)
    with pytest.raises(ValueError,match='fresh directory'):
        verify_products(tmp_path,'K00001','long',['filtered'])
