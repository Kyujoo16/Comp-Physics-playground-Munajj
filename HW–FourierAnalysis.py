import argparse
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description="Analyze and compress an astronomical light curve using Fourier Analysis."
    )
    parser.add_argument(
        "fits_file",
        type=str,
        help="Path to the input FITS file (e.g., tic0001045298.fits)",
    )
    parser.add_argument(
        "--coeffs",
        type=int,
        default=20,
        help="Number of top Fourier coefficients to keep for reconstruction (default: 20)",
    )

    args = parser.parse_args()

    hdul = fits.open("tic0001045298.fits")
    hdul.info()
    times = np.array(hdul[1].data['times'])
    fluxes = np.array(hdul[1].data['fluxes'])
    ferrs = np.array(hdul[1].data['ferrs'])

    hdul.close()


# Filtering out bad point from array
    valid_idx = ~np.isnan(times) & ~np.isnan(fluxes) & ~np.isnan(ferrs)
    times = times[valid_idx]
    fluxes = fluxes[valid_idx]
    ferrs = ferrs[valid_idx]


# Chosen continuous Epoch region
    gaps = np.diff(times)
    large_g_idx = np.argmax(gaps)

    t_s = times[:large_g_idx + 1]
    f_s = fluxes[:large_g_idx + 1]
    f_sub_detrended = f_s-np.mean(f_s)

# linear interpolation for missing data
    cadence = np.median(np.diff(t_s))
    t_uni = np.arange(t_s[0], t_s[-1], cadence)
    f_uni = np.interp(t_uni, t_s, f_sub_detrended)

# Fourier Transformation Analysis
    N = len(t_uni)
    frequencies = np.fft.rfftfreq(N, cadence)
    fft_coeff = np.fft.rfft(f_uni)

# Power spectrum Analysis
    power_spectrum = np.abs(fft_coeff) ** 2

# Chosen amount of coefficients to capture behaviour of eclipsing binary
    num_coeff_keep = 20

    sorted_indices = np.argsort(np.abs(fft_coeff))
    threshold_idx = sorted_indices[:-num_coeff_keep]

    compressed_coeff = fft_coeff.copy()
    compressed_coeff[threshold_idx] = 0

# Computation of the IFFT for signal reconsturction
    reconstructed_flux = np.fft.irfft(compressed_coeff, N)

    plt.figure(figsize=(10, 5))
    plt.errorbar(times, fluxes, fmt='b.', yerr=ferrs, markersize=2, alpha=0.5, ecolor='gray', lw=0.5, label='Raw Data')
    plt.xlabel('times(days)')
    plt.ylabel('fluxes')
    plt.title('TIC 1045298 - Raw Light Curve')
    plt.grid()
    plt.legend()

    plt.figure(figsize=(14, 10))

    plt.subplot(3, 1, 1)
    plt.plot(t_s, f_s, 'k.', alpha=0.5, label='Raw Observations')
    plt.title('Selected High-Density Epoch')
    plt.ylabel('Flux')
    plt.xlabel('Time(days)')
    plt.grid()
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(frequencies, power_spectrum, color='purple', lw=1.5)
    plt.title('Power Spectrum vs Frequency')
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.yscale('log')
    plt.grid()

    plt.subplot(3, 1, 3)
    plt.plot(t_uni, f_uni + np.mean(f_s), color='gray', alpha=0.5, label='Interpolated Grid')
    plt.plot(t_uni, reconstructed_flux + np.mean(f_s), color='blue', lw=2, label=f'Reconstructed Curve (Top {num_coeff_keep} Coeffs)')
    plt.title('Inverse Fourier Transform Reconstruction')
    plt.xlabel('Time(days)')
    plt.ylabel('Flux')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()