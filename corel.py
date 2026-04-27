import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def corel(file1):
    # 1. Read FITS file and convert to float
    with fits.open(file1) as hdul:
        img1 = hdul[0].data.astype(float)

    nx = img1.shape[1] # Note: Python is (row, col), so img1.shape[1] is width

    # 2. Remove background
    img1 -= np.median(img1)

    # 3. FFT-based Autocorrelation (Wiener-Khinchin theorem)
    ft1 = np.fft.fft2(img1)
    # conj(ft1) * ft1 is the power spectrum
    ccf = np.fft.ifft2(ft1 * np.conj(ft1)).real
    
    # Shift the zero-frequency component to the center (equivalent to IDL's shift nx/2)
    ccf = np.fft.fftshift(ccf)
    ccf -= np.median(ccf)

    # 4. Define the central zone
    m = min(256, nx // 2)
    center = nx // 2
    
    # Slice the central region: ccf1 = ccf[center-m : center+m, center-m : center+m]
    ccf1 = ccf[center-m : center+m, center-m : center+m].copy()

    # 5. Create Coordinate Grids
    # Equivalent to IDL's findgen and replicate/# 
    vals = np.arange(2 * m) - m
    x, y = np.meshgrid(vals, vals)

    # 6. Radial Masking
    r = np.sqrt(x**2 + y**2)
    rmask = 30.0
    ccf1[r < rmask] = 0.0

    # 7. Display (equivalent to tvscl)
    plt.imshow(ccf1, origin='lower', cmap='gray')
    plt.title(f"Correlation: {file1}")
    plt.show()

    # 8. Calculate Centroid of the peak
    cmax = np.max(ccf1)
    # Logic: pixels > 70% of max and in the upper half (y > 0)
    peak_indices = np.where((ccf1 > 0.7 * cmax) & (y > 0))
    
    if len(peak_indices[0]) > 0:
        tmp = ccf1[peak_indices]
        # Weighted average for centroiding
        xc = np.sum(x[peak_indices] * tmp) / np.sum(tmp)
        yc = np.sum(y[peak_indices] * tmp) / np.sum(tmp)
        
        d = np.sqrt(xc**2 + yc**2)
        print(f"Peak distance [pix]: {d:.4f}")
    else:
        print("No peak found in the upper half.")

    return d

# Example usage:
# corel('your_file.fits')

