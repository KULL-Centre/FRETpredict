""" Script for calculation of <R0> of every chromophore pair """

import numpy as np
import pandas as pd
import itertools

# Import chromophore data (Type, chromophore name, Extinction coefficient, Quantum yield)
chromophore_data = pd.read_csv('lib/R0/Dyes_extinction_QD.csv', delimiter=',', on_bad_lines='skip',
                               names=['Type', 'Chromophore', 'Ext_coeff', 'QD'])

# List of all the donors (with QD different from 0, and acceptors) and the producer
donors = chromophore_data.loc[chromophore_data['QD'] != 0, 'Chromophore'].to_list()
donors_types = chromophore_data.loc[chromophore_data['QD'] != 0, 'Type'].to_list()

print(donors, donors_types)
acceptors = chromophore_data['Chromophore'].to_list()
acceptors_types = chromophore_data['Type'].to_list()

donor_types, acceptor_types, donor, acceptor, R0 = [], [], [], [], []

donor_acceptor_types = [item for item in itertools.product(donors_types, acceptors_types)]

# Calculate R0 for every chromophore pair
for i, dye_pair in enumerate(itertools.product(donors, acceptors)):

    # Donor data (Type and emission spectrum)
    donor_type = donor_acceptor_types[i][0]

    donor_spectrum = pd.read_csv('lib/R0/{}{}.csv'.format(donor_type, dye_pair[0]))
    donor_spectrum[['Emission', 'Excitation']] = donor_spectrum[['Emission', 'Excitation']] / 100

    # Acceptor data (Type and absorption spectrum)
    acceptor_type = donor_acceptor_types[i][1]
    acceptor_spectrum = pd.read_csv('lib/R0/{}{}.csv'.format(acceptor_type, dye_pair[1]))
    acceptor_spectrum[['Emission', 'Excitation']] = acceptor_spectrum[['Emission', 'Excitation']] / 100

    # R0 calculation factors (prefactor for Angstrom, k2 for freely rotating dye, water refraction index)
    factor = 0.02108
    k2 = 2/3
    n4 = 1.4**4

    # Donor Quantum Yield
    QD = chromophore_data['QD'].loc[(chromophore_data['Chromophore'] == dye_pair[0]) &
                                     (chromophore_data['Type'] == donor_type)].item()

    # Extinction coefficient spectrum of the acceptor

    ext_coeff_max = chromophore_data['Ext_coeff'].loc[(chromophore_data['Chromophore'] == dye_pair[1]) &
                                                      (chromophore_data['Type'] == acceptor_type)].item()
    ext_coeff_acceptor = (ext_coeff_max * acceptor_spectrum['Excitation']).fillna(0)

    # Calculation of the donor spectrum integral
    donor_spectrum_integral = np.trapz(donor_spectrum['Emission'], x=donor_spectrum['Wavelength'])

    # Overlap integral between donor and acceptor
    J = np.trapz(donor_spectrum['Emission'] * ext_coeff_acceptor * donor_spectrum['Wavelength']**4,
                 x=donor_spectrum['Wavelength']) / donor_spectrum_integral

    # <R0> for the donor-acceptor dye pair
    R0_tmp = factor * np.power(k2 * QD/n4 * J, 1/6)

    if R0_tmp != 0.0:
        donor_types.append(donor_type)
        acceptor_types.append(acceptor_type)
        donor.append(dye_pair[0])
        acceptor.append(dye_pair[1])
        R0.append(R0_tmp)

    print(dye_pair, donor_type, acceptor_type, R0_tmp)

# Write data to .csv file
pd.DataFrame({'donor_type': donor_types, 'acceptor_type': acceptor_types, 'donor': donor,
              'acceptor': acceptor, 'R0': R0}).to_csv('lib/R0/R0_pairs.csv', sep='\t', index=False)

