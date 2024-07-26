import datetime
from ase import Atoms
import numpy as np
from tqdm import tqdm
from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
import pickle
import sys
from multiprocessing import Pool, Process, SimpleQueue
import os

convert_ase_to_bar=1.602176634e-19/1e-30/1e5

def generate_displacements(ph3,cutoff_pair_distance,disp_filename="phono3py_disp.yaml"):
    ph3.generate_displacements(cutoff_pair_distance=cutoff_pair_distance)
    ph3.save(disp_filename)
    print("Generated displacements")
    return ph3


def compute_forces_fc2(ph3,calculator,show_bar=True):
    # get forces for FC2
    print("Computing forces for FC2")
    print(datetime.datetime.now())
    forces = []
    for sc in tqdm(ph3.phonon_supercells_with_displacements) if show_bar else ph3.phonon_supercells_with_displacements:
        atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)
        atoms.calc = calculator
        f = atoms.get_forces()
        forces.append(f)
    
    # append 2nd order forces 
    ph3.phonon_forces = np.array(forces)
    return ph3

def compute_forces_fc3(ph3,calculator,show_bar=True):
    # get forces for FC3
    print("Computing forces for FC3")
    print(datetime.datetime.now())
    forces = []
    nat = len(ph3.supercells_with_displacements[0])
    for sc in tqdm(ph3.supercells_with_displacements) if show_bar else ph3.supercells_with_displacements:
        if sc is not None:
            atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)
            atoms.calc =calculator
            f = atoms.get_forces()
        else:
            f = np.zeros((nat, 3))
        forces.append(f)
    # append 3rd order forces
    ph3.forces = np.array(forces)
    return ph3

def _parallel_queue_worker(que, nat, calculator_function, supercell, positions_array, positions_indeces):
    calculator = calculator_function() # initialise the calculator using the provided generating function
    atoms = Atoms(supercell.symbols, cell = supercell.cell, positions=supercell.positions, pbc=True)
    atoms.calc = calculator
    for pos_, ind_ in zip(positions_array, positions_indeces):
        if np.any(np.isnan(pos_)):
            que.put([ind_, np.zeros((nat, 3))])
            continue
        atoms.set_positions(pos_)
        f = atoms.get_forces()
        que.put([ind_, f])

'''
compute the fc3 forces
inputs: ph3 Phono3py object
calculator_function: a FUNCTION which takes no arguments and generates a calculator object
no_proc: how many cpus to use (default: "auto" -- use all of them available)
show_bar: boolean to show a progress bar
'''
def compute_forces_fc3_parallel2(ph3, calculator_function, no_proc='auto', show_bar=True):
    print("Computing forces for FC3 (parallel v2)")
    if no_proc == 'auto':
        no_proc = int(os.getenv("OMP_NUM_THREADS", "1"))
        print(f"detected {no_proc} cpus available")
        if no_proc > 1:
            no_proc -= 1 # minus one (?) for master process
    if not isinstance(no_proc, int):
        raise ValueError()
    print(datetime.datetime.now())
    ### setup stuff
    nat = len(ph3.supercells_with_displacements[0])
    n_displacements = len(ph3.supercells_with_displacements)
    supercell0 = ph3.supercells_with_displacements[0]
    que = SimpleQueue()
    #print(f"f{nat=} {n_displacements=}")
    #### run no_proc processes with equal portions of the job
    positions_atoms = [sc.positions if not sc is None else np.full((nat, 3), np.nan) for sc in ph3.supercells_with_displacements]
    split_positions = np.array_split(positions_atoms, no_proc)
    split_indeces = np.array_split(list(range(n_displacements)), no_proc)
    processes = [Process(target=_parallel_queue_worker, args=(que, nat, calculator_function, supercell0, positions_, indeces_)) for positions_, indeces_ in zip(split_positions, split_indeces)]
    for p in processes:
        p.start()
    forces = np.zeros((n_displacements, nat, 3))
    # collect fc3 results as they come; display progress with tqdm
    for i in tqdm(range(n_displacements)) if show_bar else range(n_displacements):
        ind, f = que.get()
        forces[ind, :, :] = f
    ph3.forces = forces
    print("parallel force computation complete @ ", datetime.datetime.now(), flush=True)
    return ph3

def save_ph3(ph3,forces_output_file):
    print("Writing FC3 FORCES SET")
    print(datetime.datetime.now())

    ph3.save(forces_output_file,
                settings={'force_sets': True,
                    'displacements': True,
                    'force_constants': False,
                    'born_effective_charge': False,
                    'dielectric_constant': False})

    print("Finished")
    print(datetime.datetime.now())

def save_ph3_pickle(ph3,pickle_fname):
    with open(pickle_fname,"wb") as f:
        pickle.dump(ph3.phonon_forces,f )
        pickle.dump(ph3.forces,f)


def run_forces(ph3,calculator, cutoff_pair_distance,disp_filename="phono3py_disp.yaml",show_bar=True):
    generate_displacements(ph3,cutoff_pair_distance,disp_filename)
    compute_forces_fc2(ph3,calculator,show_bar=show_bar)
    compute_forces_fc3(ph3,calculator,show_bar=show_bar)
    return ph3

def run_fc2_forces(ph3,calculator, cutoff_pair_distance,disp_filename="phono3py_disp.yaml",show_bar=True):
    generate_displacements(ph3,cutoff_pair_distance,disp_filename)
    compute_forces_fc2(ph3,calculator,show_bar=show_bar)
    write_ph3_fc2_hdf5(ph3)
    return ph3

def write_ph3_fc2_hdf5(ph3):
    ph3.produce_fc2()
    print("Writing FC2 HDF5")
    print(datetime.datetime.now())
    write_fc2_to_hdf5(
        ph3.fc2,
        p2s_map=ph3.phonon_primitive.p2s_map,
        physical_unit="eV/angstrom^2",
    )
    print("Finished writing FC2 HDF5")
    print(datetime.datetime.now())

def write_ph3_fc3_hdf5(ph3):
    print("producing fc3... @ ", datetime.datetime.now(), flush=True)
    ph3.produce_fc3()
    print("Writing FC3 HDF5")
    print(datetime.datetime.now())
    write_fc3_to_hdf5(
        ph3.fc3,
        p2s_map=ph3.phonon_primitive.p2s_map,
    )
    print("Finish writing FC3 HDF5")
    print(datetime.datetime.now(), flush=True)

def print_stresses(ph3,calculator,output=sys.stdout):
    sc=ph3.unitcell
    atoms = Atoms(sc.symbols, cell=sc.cell, positions=sc.positions, pbc=True)
    atoms.calc=calculator
    print("Final Energy", atoms.get_potential_energy()," ev",file=output)
    print("Final Stress",atoms.get_stress()*convert_ase_to_bar," bar",file=output)

