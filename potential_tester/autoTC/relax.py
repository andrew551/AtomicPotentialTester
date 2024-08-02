import numpy as np
from ase.io import read, write
from ase.io.vasp import read_vasp
from ase import Atoms
from ase.calculators.lammpslib import LAMMPSlib
from ase.constraints import FixSymmetry
from ase.filters import UnitCellFilter, ExpCellFilter, StrainFilter,FrechetCellFilter
from ase.optimize import BFGS, FIRE, MDMin, GPMin
from ase.spacegroup.symmetrize import check_symmetry
from potential_tester.load import *

import warnings, os




convert_ase_to_bar=1.602176634e-19/1e-30/1e5
fmax_thresholds_default=np.array([1e1,5e-1,1e-1,5e-2,1e-2,5e-3,1e-3,1e-3,1e-3,5e-4,5e-4,5e-4,1e-4,1e-4,5e-5,5e-5,1e-5,1e-5,1e-5,1e-5,8e-6,8e-6,5e-6,5e-6,3e-6,3e-6])
Optimizer_default=BFGS



def vc_fc_relax(atoms,num_iter,thresholds,Optimizer=Optimizer_default,allow_tilt=True,alpha=None):
    vc_mask=None
    if not allow_tilt:
        vc_mask=[True,True,True,False,False,False]


    #fc_filter=ExpCellFilter(atoms, mask=[False]*6)
    #vc_filter=StrainFilter(atoms,mask=vc_mask)
    for i in range(num_iter):
        #vc relax
        vc_filter=StrainFilter(atoms,mask=vc_mask)
        if Optimizer.__name__ in ["BFGS",'SciPyFminCG']:
            dyn = Optimizer(vc_filter,alpha=alpha)
        else:
            dyn = Optimizer(vc_filter)
        dyn.run(fmax=thresholds[i])
        print("Energy", atoms.get_potential_energy()," ev")
        print("Stress",atoms.get_stress()*convert_ase_to_bar," bar")
        
        #fc-relax (fixed cell)
        fc_filter=ExpCellFilter(atoms, mask=[False]*6)
        dyn = Optimizer(fc_filter)
        dyn.run(fmax=thresholds[i])
        print("Energy", atoms.get_potential_energy()," ev")
        print("Stress",atoms.get_stress()*convert_ase_to_bar," bar")
    return atoms

def reset_threshold(atoms,thresholds):
    fmax=np.sqrt((atoms.get_forces() ** 2).sum(axis=1).max())
    return thresholds[fmax>thresholds]

def run_rvcrelax(atoms,calculator,fmax_thresholds=fmax_thresholds_default,threshold_min=1e-10,allow_tilt=False,Optimizer=Optimizer_default,alpha=None):
    print("Initial symmetry at precision 1e-6")
    check_symmetry(atoms, 1.0e-6, verbose=True)

    fmax_thresholds=np.array(fmax_thresholds).copy()
    fmax_thresholds=fmax_thresholds[fmax_thresholds>threshold_min]

    # Now we repeat the optimisation with the symmetrization constraint in place
    atoms.set_constraint(FixSymmetry(atoms))
    atoms.calc=calculator


    print('Reseting thresholds to be lower than maximal force: ')
    fmax_thresholds=reset_threshold(atoms,fmax_thresholds)
    num_iter=len(fmax_thresholds)
    print(fmax_thresholds)

    print("Initial Energy", atoms.get_potential_energy()," ev")
    print("Initial Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    if len(fmax_thresholds)==0:
        fmax_thresholds=fmax_thresholds_default[fmax_thresholds_default>threshold_min]

    vc_fc_relax(atoms,num_iter,fmax_thresholds,alpha=alpha)

    print("After keeping symmetry VC/FC relax Energy", atoms.get_potential_energy()," ev")
    print("After keeping symmetry VC/FC relax Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    # We print out the initial symmetry groups at two different precision levels
    print("After keeping symmetry VC/FC relax symmetry at precision 1e-5")
    check_symmetry(atoms, 1.0e-5, verbose=True)

    print("Minimising without keeping symmetries")

    del atoms.constraints
    if len(fmax_thresholds)==0:
        fmax_thresholds=[1e-5]
        print("Thresholds empty! Reset to 1e-5.")
    fmax_thresholds=[fmax_thresholds[-1]]*4
    num_iter=len(fmax_thresholds)
    vc_fc_relax(atoms,num_iter,fmax_thresholds,allow_tilt=allow_tilt)
    



    print("Final Energy", atoms.get_potential_energy()," ev")
    print("Final Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    print("Final symmetry at precision 1e-4")
    check_symmetry(atoms, 1.0e-4, verbose=True)
    print("Final symmetry at precision 1e-5")
    check_symmetry(atoms, 1.0e-5, verbose=True)
    

    atoms_write=Atoms(symbols=atoms.symbols,positions=atoms.positions,cell=atoms.cell,pbc=atoms.pbc)
    write("POSCAR",atoms_write,format='vasp')
    return atoms

def run_joint_relax(atoms,calculator,fmax=1e-4,fmax_final=None,allow_tilt=False,Optimizer=Optimizer_default,alpha=None,force_symmetry=True):
    if fmax_final is None:
        fmax_final=fmax


    atoms.calc=calculator
    
    no_tilt_mask=[True,True,True,False,False,False]
    vc_mask=None
    if not allow_tilt:
        vc_mask=no_tilt_mask
    
    print("Initial Energy", atoms.get_potential_energy()," ev")
    print("Initial Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    print("Initial symmetry at precision 1e-6")
    sym_before_5=check_symmetry(atoms, 1.0e-5, verbose=True)
    sym_before_3=check_symmetry(atoms, 1.0e-3, verbose=True)

    atoms.set_constraint(FixSymmetry(atoms))

    vc_filter=StrainFilter(atoms,mask=vc_mask)
    exp_filter=ExpCellFilter(atoms,mask=no_tilt_mask)
    if Optimizer.__name__ in ["BFGS",'SciPyFminCG']:
        dyn_cell = Optimizer(vc_filter,alpha=alpha)
        dyn_atoms_only = Optimizer(atoms,alpha=alpha)
        dyn_total=Optimizer(exp_filter,alpha=alpha)
    else:
        dyn_cell = Optimizer(vc_filter)
        dyn_atoms_only = Optimizer(atoms)
        dyn_total=Optimizer(exp_filter)
    


    # Run a optimisation for atomic positions 
    # with every step rescaling the cell to minimise stress
    #dyn_total.run(fmax=1e-3,steps=200)
    
    for _ in dyn_atoms_only.irun(fmax=fmax,steps=500):
        dyn_cell.run(fmax=fmax,steps=1000)
    dyn_atoms_only.run(fmax=fmax_final,steps=500)

    print("After keeping symmetry VC/FC relax Energy", atoms.get_potential_energy()," ev")
    print("After keeping symmetry VC/FC relax Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    # We print out the initial symmetry groups at two different precision levels
    print("After keeping symmetry VC/FC relax symmetry at precision 1e-5")
    sym_middle_5=check_symmetry(atoms, 1.0e-5, verbose=True)

    if sym_middle_5['number']!=sym_before_5['number']:
        warnings.warn(f"SYMMETRY IS NOT KEPT DURING FxSymmetry RELAXTION in folder {os.getcwd()}")

    # delete constrainsts and run a optimisation for atomic positions 
    # with every step rescaling the cell to minimise stress
    atoms_symmetry=atoms.copy()

    del atoms.constraints

    dyn_atoms_only.run(fmax=fmax_final,steps=200)

    print("Right after deleting symmetry VC/FC relax Energy", atoms.get_potential_energy()," ev")
    print("Right after deleting symmetry VC/FC relax Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    # We print out the initial symmetry groups at two different precision levels
    print("Right after deleting symmetry VC/FC relax symmetry at precision 1e-5")
    check_symmetry(atoms, 1.0e-5, verbose=True)
    check_symmetry(atoms, 1.0e-3, verbose=True)


    for _ in dyn_atoms_only.irun(fmax=fmax,steps=200):
        dyn_cell.run(fmax=fmax,steps=500)
    dyn_atoms_only.run(fmax=fmax_final,steps=200)

    print("Final Energy", atoms.get_potential_energy()," ev")
    print("Final Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    print("Final symmetry at precision 1e-4")
    check_symmetry(atoms, 1.0e-4, verbose=True)
    print("Final symmetry at precision 1e-5")
    sym_after_5=check_symmetry(atoms, 1.0e-5, verbose=True)
    

    # compare symmetries
    
    if sym_middle_5['number']!=sym_after_5['number'] and force_symmetry:
        atoms=atoms_symmetry
        warnings.warn(f"SYMMETRY IS NOT KEPT AFTER DELETING CONSTRAINT, redirecting to structure with symmetry, in folder {os.getcwd()}")

    atoms_write=Atoms(symbols=atoms.symbols,positions=atoms.positions,cell=atoms.cell,pbc=atoms.pbc)
    write("POSCAR",atoms_write,format='vasp')
    return atoms

    

def run_novc_relax(atoms,calculator,fmax=1e-4,fmax_final=None,allow_tilt=False,Optimizer=Optimizer_default,alpha=None,force_symmetry=True):
    if fmax_final is None:
        fmax_final=fmax


    atoms.calc=calculator
    
    
    print("Initial Energy", atoms.get_potential_energy()," ev")
    print("Initial Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    print("Initial symmetry at precision 1e-6")
    sym_before_5=check_symmetry(atoms, 1.0e-5, verbose=True)
    sym_before_3=check_symmetry(atoms, 1.0e-3, verbose=True)

    atoms.set_constraint(FixSymmetry(atoms))

    if Optimizer.__name__ in ["BFGS",'SciPyFminCG']:
        dyn_atoms_only = Optimizer(atoms,alpha=alpha)
    else:
        dyn_atoms_only = Optimizer(atoms)
    


    # Run a optimisation for atomic positions 
    # with every step rescaling the cell to minimise stress
    #dyn_total.run(fmax=1e-3,steps=200)
    
    dyn_atoms_only.run(fmax=fmax_final,steps=500)

    print("After keeping symmetry VC/FC relax Energy", atoms.get_potential_energy()," ev")
    print("After keeping symmetry VC/FC relax Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    # We print out the initial symmetry groups at two different precision levels
    print("After keeping symmetry VC/FC relax symmetry at precision 1e-5")
    sym_middle_5=check_symmetry(atoms, 1.0e-5, verbose=True)

    if sym_middle_5['number']!=sym_before_5['number']:
        warnings.warn(f"SYMMETRY IS NOT KEPT DURING FxSymmetry RELAXTION in folder {os.getcwd()}")

    # delete constrainsts and run a optimisation for atomic positions 
    # with every step rescaling the cell to minimise stress
    atoms_symmetry=atoms.copy()

    del atoms.constraints

    dyn_atoms_only.run(fmax=fmax_final,steps=200)

    print("Final Energy", atoms.get_potential_energy()," ev")
    print("Final Stress",atoms.get_stress()*convert_ase_to_bar," bar")

    print("Final symmetry at precision 1e-4")
    check_symmetry(atoms, 1.0e-4, verbose=True)
    print("Final symmetry at precision 1e-5")
    sym_after_5=check_symmetry(atoms, 1.0e-5, verbose=True)
    

    # compare symmetries
    
    if sym_middle_5['number']!=sym_after_5['number'] and force_symmetry:
        atoms=atoms_symmetry
        warnings.warn(f"SYMMETRY IS NOT KEPT AFTER DELETING CONSTRAINT, redirecting to structure with symmetry, in folder {os.getcwd()}")

    atoms_write=Atoms(symbols=atoms.symbols,positions=atoms.positions,cell=atoms.cell,pbc=atoms.pbc)
    write("POSCAR",atoms_write,format='vasp')
    return atoms



def run_rvcrelax_ph3(ph3,calculator,fmax_thresholds=fmax_thresholds_default,threshold_min=1e-10,allow_tilt=False,Optimizer=Optimizer_default,alpha=None):
    atoms=phono3py2ase(ph3)
    atoms=run_rvcrelax(atoms,calculator,fmax_thresholds=fmax_thresholds_default,threshold_min=threshold_min,allow_tilt=allow_tilt,Optimizer=Optimizer,alpha=alpha)
    ph3=phono3py_ase2phono3py(ph3,atoms)
    return ph3

def run_relax_ph3(ph3,calculator,fmax=1e-4,fmax_final=None,fmax_thresholds=fmax_thresholds_default,threshold_min=1e-10,allow_tilt=False,Optimizer=Optimizer_default,alpha=None,force_symmetry=True):
    atoms=phono3py2ase(ph3)
    atoms=run_joint_relax(atoms,calculator,fmax=fmax,fmax_final=fmax_final,allow_tilt=allow_tilt,Optimizer=Optimizer,alpha=alpha,force_symmetry=force_symmetry)
    ph3=phono3py_ase2phono3py(ph3,atoms)
    return ph3

def run_CHGNet_relax_ph3(ph3,output_POSCAR="POSCAR"):
    from autoTC.CHGNet_driver import relax_CHGNet_ASE
    atoms=phono3py2ase(ph3)
    atoms=relax_CHGNet_ASE(atoms)
    ph3=phono3py_ase2phono3py(ph3,atoms)
    write(output_POSCAR,atoms,format='vasp')
    return ph3

def run_M3GNet_relax_ph3(ph3,output_POSCAR="POSCAR"):
    from autoTC.CHGNet_driver import relax_M3GNet_ASE
    atoms=phono3py2ase(ph3)
    atoms=relax_M3GNet_ASE(atoms)
    write(output_POSCAR,atoms,format='vasp')
    ph3=phono3py_ase2phono3py(ph3,atoms)
    return ph3