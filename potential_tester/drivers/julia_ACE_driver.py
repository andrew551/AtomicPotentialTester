from ase.calculators.lammpslib import LAMMPSlib

def setup_ACE_lammps_calculator_SiO(path_ACE_potential):

    lammps_cmds = ['pair_style      hybrid/overlay pace table linear 6000',
        f'pair_coeff      * * pace {path_ACE_potential}.yace O Si',
        f'pair_coeff      1 1 table {path_ACE_potential}_pairpot.table O_O',
        f'pair_coeff      1 2 table {path_ACE_potential}_pairpot.table O_Si',
        f'pair_coeff      2 2 table {path_ACE_potential}_pairpot.table Si_Si' 
    ]
    
    lammps_header=[ 'units metal',
                    'boundary p p p',
                    'box tilt large',
                    'atom_style atomic',
                    'atom_modify map array sort 0 0'
        ]

    atoms_dict={'O':1,'Si':2}

    lammps_log_file=None#'lammps_out.log'


    #lammps calculator
    calculator=LAMMPSlib(lmpcmds=lammps_cmds,
                        log_file=lammps_log_file,
                        lammps_header=lammps_header, 
                        atom_types=atoms_dict,
                        keep_alive=True)

    return calculator
