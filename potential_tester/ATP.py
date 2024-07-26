import sys
import os
from pathlib import Path
from potential_tester.drivers import julia_ACE_driver 
from potential_tester import file_conversion
from potential_tester import load
from potential_tester.autoTC import fc_interface, relax, tc_calc
import phonopy
import matplotlib.pyplot as plt
import numpy as np
import h5py

class ATP:

    def __init__(self, **kwaargs):
        self.config = {
            'potential_type': 'ACE_julia',
            'potential_path' : '',
            'workdir' : os.getcwd(),
            'atoms_file': 'POSCAR',
            'supercell': [4,4,4],
            'FC2_SC':[4,4,4],
            'FC3_SC':[4,4,4],
            'Q_MESH':[11,11,11],
            'save_graphics':False,
            'save_force_constants':True,
            'show_progress_bar':True,
        }
        self.update_config_parameters(kwaargs)
        self.initialise_calculator()
        self.atoms = file_conversion.read_reg(self.config['atoms_file'])
        self.ph3 = load.ase2phono3py(self.atoms, self.config['FC2_SC'], self.config['FC3_SC'])



    def update_config_parameters(self, parameters):
        self.config.update(parameters)

    def initialise_calculator(self):
        if self.config['potential_type'] == 'ACE_julia':
            self.calculator = julia_ACE_driver.setup_ACE_lammps_calculator_SiO(self.config['potential_path'])
            
        else:
            raise Exception("unimplemented")

    def run_relax(self):
        initial_density = file_conversion.get_density(self.atoms)
        self.atoms = relax.run_joint_relax(self.atoms,self.calculator)
        self.ph3 = load.ase2phono3py(self.atoms, self.config['FC2_SC'], self.config['FC3_SC'])
        file_conversion.write_file(self.atoms, Path(self.config['workdir']) / 'POSCAR_relaxed', out_type='vasp')
        final_density = file_conversion.get_density(self.atoms)
        with open(Path(self.config['workdir']) / 'density_results.txt', 'w') as f:
            f.writelines(['initial density: ' + str(initial_density) + '\n', 'final density:' + str(final_density) + '\n'])

    # compute fc2; save data
    def compute_force_constants(self):
        fc_interface.run_forces(self.ph3,self.calculator, cutoff_pair_distance=None,disp_filename="phono3py_disp.yaml")
        fc_interface.write_ph3_fc2_hdf5(self.ph3)
        fc_interface.write_ph3_fc3_hdf5(self.ph3)
        pickle_filename="forces.pckl"
        fc_interface.save_ph3_pickle(self.ph3, pickle_filename)
    
    # load ph3 fc2/fc3 if exists
    def load_ph3(self):
        #raise Exception("unimplemented")
        workdir = Path(self.config['workdir']) 
        input_file=workdir/"phono3py_disp.yaml"
        ph3=load.load_ph3py_yaml(input_file)

        pickle_filename=workdir/"forces.pckl"
        ph3=load.load_pickled_forces(ph3,pickle_filename)
        self.atoms = load.phono3py2ase(ph3)
    
    def compute_specific_heat(self):
        phonon = self.as_phonopy()
        phonon.run_mesh(self.config['Q_MESH'])
        phonon.run_thermal_properties(t_step=10,
                                    t_max=1000,
                                    t_min=0)
        tp_dict = phonon.get_thermal_properties_dict()
        temperatures = tp_dict['temperatures']
        free_energy = tp_dict['free_energy']
        entropy = tp_dict['entropy']
        heat_capacity = tp_dict['heat_capacity']
        fig, ax = plt.subplots()
        plt.plot(temperatures, heat_capacity)
        plt.savefig(Path(self.config['workdir']) / 'heat_capacity.png')


    def compute_phonons(self):
        phonon = self.as_phonopy()
        phonon.auto_band_structure(plot=True)
        plt.savefig(Path(self.config['workdir']) / 'phonons.png')

    def compute_tc(self, T_min=80, T_max=640, T_step=80):
        tc_calc.calc_tc_hdf5_yaml(T_min,T_max,T_step,initial_string="phono3py_disp.yaml",output_file="tc.out",Q_MESH=[9,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py="/mnt/userapps/q13camb_apps/python/packages/phono3py-OMP")
        # plot tc
        plot_tc_data(file_in=Path(self.config['workdir']) / 'kappa-m999.hdf5')

    def compute_gruneisen(self):
        raise Exception("unimplemented")
    
    def as_phonopy(self):
        # return phonopy version of self
        phonon = phonopy.load(supercell_matrix=self.config['FC2_SC'], unitcell_filename='POSCAR_relaxed', force_constants_filename='fc2.hdf5')
        return phonon
        
    def compute_bulk_modulus(self):
        # need fc-relax
        raise Exception("unimplemented")

def plot_tc_data(file_in = 'kappa-m999.hdf5'):
    f = h5py.File(file_in,'r')

    kappa_tot = np.mean(f['kappa_TOT_RTA'][:, :3], axis=1) # (kappa_xx + kappa_yy + kappa_zz) / 3
    kappa_P = np.mean(f['kappa_P_RTA'][:, :3], axis=1) # (kappa_xx + kappa_yy + kappa_zz) / 3
    kappa_C = np.mean(f['kappa_C'][:, :3], axis=1) # (kappa_xx + kappa_yy + kappa_zz) / 3
    temperatures = np.array(f['temperature'][:])


    labels = ['$\\kappa_{tot}$', '$\\Gamma_{rel}=0.5$', '$\\Gamma_{rel}=2$', '$\\kappa_P$', '$\\kappa_C$']

    fig, ax = plt.subplots()

    #ax.fill_between(temperatures, k_half, k_double, color='lightgray')
    ax.plot(temperatures, kappa_P, color='#FAA43A', label = labels[3], linewidth=2)
    ax.plot(temperatures, kappa_C, color='#60BD68', label = labels[4], linewidth=2)
    ax.plot(temperatures, kappa_tot, color='black', label = labels[0], linewidth=3)

    #plt.ylim(0, 2)
    #plt.xlim(60, 440)
    plt.legend(loc=(0.79, 0.25), frameon=False, fontsize=14)
    plt.tick_params(direction='in')
    plt.xlabel('$T$ (K)', fontsize=15)
    plt.ylabel('$\\kappa(T)$ $\\mathrm{(W m^{-1} K^{-1})}$', fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    #plt.yticks(np.linspace(0, 2, 5))
    plt.savefig("conductivity_plot.png", dpi=600, bbox_inches='tight')