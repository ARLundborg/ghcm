import scipy.io
import pandas as pd
import numpy as np

mat = scipy.io.loadmat("../extdata/original.mat")

y_df = pd.DataFrame(mat["y"], columns=["date", "color", "ash"]).drop(axis=1, labels="date")

y_df.to_csv("../../data-raw/ash_color.txt", header=True, index=False)

emission_spectra = [mat["X"][:, i:(i+571)] for i in range(7)]
emission_wavelengths = np.linspace(275, 560, 571)

excitation_340_df = pd.DataFrame(emission_spectra[0], columns=emission_wavelengths)
excitation_340_df.to_csv("excitation_340.txt", header=True, index=False, float_format='%.3f')

excitation_325_df = pd.DataFrame(emission_spectra[1], columns=emission_wavelengths)
excitation_325_df.to_csv("excitation_325.txt", header=True, index=False, float_format='%.3f')

excitation_305_df = pd.DataFrame(emission_spectra[2], columns=emission_wavelengths)
excitation_305_df.to_csv("excitation_305.txt", header=True, index=False, float_format='%.3f')

excitation_290_df = pd.DataFrame(emission_spectra[3], columns=emission_wavelengths)
excitation_290_df.to_csv("excitation_290.txt", header=True, index=False, float_format='%.3f')

excitation_255_df = pd.DataFrame(emission_spectra[4], columns=emission_wavelengths)
excitation_255_df.to_csv("excitation_255.txt", header=True, index=False, float_format='%.3f')

excitation_240_df = pd.DataFrame(emission_spectra[5], columns=emission_wavelengths)
excitation_240_df.to_csv("excitation_240.txt", header=True, index=False, float_format='%.3f')

excitation_230_df = pd.DataFrame(emission_spectra[6], columns=emission_wavelengths)
excitation_230_df.to_csv("excitation_230.txt", header=True, index=False, float_format='%.3f')
