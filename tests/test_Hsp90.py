import pytest
import MDAnalysis
import numpy as np
from FRETpredict.FRET import FRETpredict

def test_HIV1PR_PDBs():
    u = MDAnalysis.Universe('tests/data/HIV-1PR/PDBs/HIV-1PR.pdb')
    DEER = DEERpredict(u, residues = [55, 55], chains=['A', 'B'], log_file = 'tests/data/HIV-1PR/log', rmax=7.5)
    DEER.run(output_prefix = 'tests/data/HIV-1PR/res', filter_stdev = 0.0707)
    r, p = np.loadtxt('tests/data/HIV-1PR/res-55-55.dat',unpack=True)
    r_ref, p_ref = np.loadtxt('tests/data/HIV-1PR/DEER_HIV-1PR.dat',unpack=True)
    p_ref /= np.trapz(p_ref, r_ref)
    assert np.power(p-p_ref,2).sum() < 0.001
    for pdb in ['3bvb','2bpx','1hhp','1tw7']:
        u = MDAnalysis.Universe('tests/data/HIV-1PR/PDBs/{:s}.pdb'.format(pdb))
        DEER = DEERpredict(u, residues = [55,55], chains=['A', 'B'], rmax=7.5)
        DEER.run(output_prefix = 'tests/data/HIV-1PR/PDBs/'+pdb, filter_stdev = 0.0707, dt = 5.49/511)
        x, y = np.loadtxt('tests/data/HIV-1PR/PDBs/'+pdb+'-55-55_time-domain.dat',unpack=True)
        x_ref, y_ref = np.loadtxt('tests/data/HIV-1PR/PDBs/REF'+pdb+'-55-55_time-domain.dat',unpack=True)
        assert np.power(y-y_ref,2).sum() < 0.0001
        r, p = np.loadtxt('tests/data/HIV-1PR/PDBs/'+pdb+'-55-55.dat',unpack=True)
        r_ref, p_ref = np.loadtxt('tests/data/HIV-1PR/PDBs/REF'+pdb+'-55-55.dat',unpack=True)
        assert np.power(p-p_ref,2).sum() < 0.0001

def test_HIV1PR_SIMs():
    for sim in ['unbiased','rdc']:
        u = MDAnalysis.Universe('tests/data/HIV-1PR/sims/{:s}.pdb'.format(sim),'tests/data/HIV-1PR/sims/{:s}.xtc'.format(sim))
        DEER = DEERpredict(u, residues = [55,55], chains=['A', 'B'], z_cutoff = 0.05, rmax=7.5)
        DEER.run(output_prefix='tests/data/HIV-1PR/sims/{:s}'.format(sim), filter_stdev = 0.0707, dt = 5.49/511)
        x, y = np.loadtxt('tests/data/HIV-1PR/sims/'+sim+'-55-55_time-domain.dat',unpack=True)
        x_ref, y_ref = np.loadtxt('tests/data/HIV-1PR/sims/REF'+sim+'-55-55_time-domain.dat',unpack=True)
        assert np.power(y-y_ref,2).sum() < 0.0001
        r, p = np.loadtxt('tests/data/HIV-1PR/sims/'+sim+'-55-55.dat',unpack=True)
        r_ref, p_ref = np.loadtxt('tests/data/HIV-1PR/sims/REF'+sim+'-55-55.dat',unpack=True)
        assert np.power(p-p_ref,2).sum() < 0.0001

def test_T4L():
    for pdb in ['3dmv','2lcb','2lc9']:
        u = MDAnalysis.Universe('tests/data/T4L/PDBs/{:s}.pdb'.format(pdb))
        for residues in [[89, 109],[109, 140]]:
            DEER = DEERpredict(u, residues = residues, temperature = 298, z_cutoff = 0.05, rmax = 7.5)
            DEER.run(output_prefix = 'tests/data/T4L/PDBs/{:s}'.format(pdb), filter_stdev = 0.0707)
            r, p = np.loadtxt('tests/data/T4L/PDBs/{:s}-{:d}-{:d}.dat'.format(pdb,*residues),unpack=True)
            r_ref, p_ref = np.loadtxt('tests/data/T4L/PDBs/REF{:s}-{:d}-{:d}.dat'.format(pdb,*residues),unpack=True)
            assert np.power(p-p_ref,2).sum() < 0.0001
