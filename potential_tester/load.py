from phono3py.cui.load import load
from ase.io import read, write
from ase.io.vasp import read_vasp
from ase import Atoms
from ase.spacegroup.symmetrize import check_symmetry
from phonopy.structure.atoms import PhonopyAtoms
from phono3py.api_phono3py import Phono3py
from phonopy.api_phonopy import Phonopy
from phonopy import file_IO
import h5py

import pickle

def load_ph3py_yaml(input_file,produce_fc=False,symmetrize_fc=False):
    ph3=load(input_file,produce_fc=produce_fc,symmetrize_fc=symmetrize_fc)
    return ph3


def load_pickled_forces(ph3,pickle_filename):
    with open(pickle_filename,"rb") as f:
        ph3.phonon_forces=pickle.load(f)
        ph3.forces=pickle.load(f)
    
    return ph3


def phono3py2ase(ph3):
    phonoatoms=ph3.unitcell
    atoms = Atoms(phonoatoms.symbols, cell=phonoatoms.cell, positions=phonoatoms.positions, pbc=True)
    return atoms

def ase2phonoatoms(atoms):
    phonoatoms = PhonopyAtoms(atoms.symbols, cell=atoms.cell, positions=atoms.positions, pbc=True)
    return phonoatoms

def phono3py_ase2phono3py(ph3,atoms):
    unitcell = ase2phonoatoms(atoms)
    pscm=ph3.phonon_supercell_matrix
    scm=ph3.supercell_matrix
    pm=ph3.primitive_matrix
    ph3_new=Phono3py(unitcell=unitcell,supercell_matrix=scm,phonon_supercell_matrix=pscm,primitive_matrix=pm)
    #ph3_new.nac_params(ph3.nac_params)

    return ph3_new

def phono3py2phonopy(ph3):
    return Phonopy(unitcell=ph3.unitcell,primitive_matrix=ph3.primitive_matrix,supercell_matrix=ph3.phonon_supercell_matrix)

def phonopy_load_fc2(phonopy,filename="fc2.hdf5"):
    force_constants    = file_IO.read_force_constants_hdf5( filename=filename )
    phonopy.set_force_constants(force_constants)
    phonopy.symmetrize_force_constants()  
    return phonopy

def write_gv_volume_hdf5(filename="gv_vol.hdf5",gv=None,volume=None):
    with h5py.File(filename,"w") as f:
        f.create_dataset("group_velocity",data=gv,compression=4)
        f.create_dataset("volume",data=volume)

    
def write_gv_volume_density_hdf5(filename="gv_vol.hdf5",gv=None,volume=None,density=None):
    with h5py.File(filename,"w") as f:
        f.create_dataset("group_velocity",data=gv,compression=4)
        f.create_dataset("volume",data=volume)
        f.create_dataset("density",data=volume)

def ase2phono3py(atoms,FC2_SC,FC3_SC):
    unitcell=ase2phonoatoms(atoms)
    return Phono3py(unitcell=unitcell,supercell_matrix=FC3_SC,phonon_supercell_matrix=FC2_SC)
