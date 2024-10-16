from ase.calculators.lammpslib import LAMMPSlib

def setup_Deringher_ACE_lammps_calculator_SiO(path_ACE_potential):

    lammps_cmds = ['pair_style  pace',
        f'pair_coeff      * * {path_ACE_potential}/SiOx_potential.yace O Si',
    ]
    
    lammps_header=[ 'units metal',
                    'boundary p p p',
                    'box tilt large',
                    'atom_style atomic',
                    'atom_modify map array sort 0 0'
        ]

    atoms_dict={'O':1,'Si':2}

    lammps_log_file='lammps_out.log'


    #lammps calculator
    calculator=LAMMPSlib(lmpcmds=lammps_cmds,
                        log_file=lammps_log_file,
                        lammps_header=lammps_header, 
                        atom_types=atoms_dict,
                        keep_alive=True)

    return calculator