import numpy as np
import matplotlib.pyplot as plt 
from scipy.misc import derivative
from pathlib import Path 

#from  potential_tester.ACE_mixer.CHIK import compute_CHIK
from CHIK import compute_CHIK
from ZBL import make_ZBL_potential
def smooth_transition_function_normalised(x, m = 2):
    x = np.array(x)
    ret = np.zeros(x.shape)
    inrange = np.logical_and(-1 < x, x < 1)
    xx = x[inrange]
    ret[inrange] = np.tanh(m * xx / (1-xx*xx))
    ret[x <= -1] = -1
    ret[x >= 1] = 1
    return ret

# m: gradient of normalised function at x=0 (must be greater than sqrt(3)): default 2
def make_transition_function(start, end, m=2):
    def stf(x):
        u = (x-start) / (end - start)
        return (smooth_transition_function_normalised(2*(u-0.5), m=m)+1) / 2
    return stf

def load_ACE_pairpot(filename):
    state = 'init'
    all_tables = {}
    with open(filename) as f:
        for line in f:
            vals = line.split()
            if not vals:
                continue # empty line
            if vals[0][0] == '#':
                continue
            if state == 'init':
                specie = vals[0]
                state = 'init2'
                table = []
            elif state == 'init2':
                if not vals[0] == 'N':
                    raise Exception("expected N")
                N = int(vals[1])
                state = 'read'
                n_read = 0
            elif state == 'read':
                table.append(list(map(float, vals)))
                n_read +=1
                if n_read == N:
                    all_tables[specie] = np.array(table)
                    state = 'init'
    return all_tables

def write_ACE_pairpot(filename, ace_table):
    with open(filename, 'w') as f:
        f.write("# DATE: none UNITS: metal CONTRIBUTOR: ACE1.jl - https://github.com/ACEsuit/ACE1.jl\n")
        f.write("# ACE1 pair potential\n\n")

        for specie, table in ace_table.items():
            f.write(specie+'\n')
            f.write(f'N {table.shape[0]}\n\n')
            for i in range(table.shape[0]):
                f.write(f'{i+1} {table[i, 1]} {table[i, 2]} {table[i, 3]}\n')
            f.write('\n')

atomic_numbers = {'Si' : 14, 'O': 8, 'C':6, 'H':1}
def created_zbl_mixed_ACE(ACE_table):
    new_ACE_table = {}
    transition_function = make_transition_function(0.7, 1.0, m=2)
    for specie, table in ACE_table.items():
        s1, s2 = specie.split('_')
        z1, z2 = atomic_numbers[s1], atomic_numbers[s2]
        zbl_potential = make_ZBL_potential(z1, z2, epsilon = 0.2)
        r = table[:, 1]
        v_zbl = zbl_potential(r)
        f_zbl = -derivative(zbl_potential, r, dx=1e-4) # F  = -dE/dr
        weight = transition_function(r)
        weight_deri = derivative(transition_function, r, dx = 1e-4)
        mixed_v = table[:, 2] * weight + v_zbl * (1 - weight)
        mixed_f = table[:, 3] * weight + f_zbl * (1 - weight) - weight_deri * (table[:, 2] - v_zbl) # remember product rule!

        new_table = np.copy(table)
        new_table[:, 2] = mixed_v
        new_table[:, 3] = mixed_f

        new_ACE_table[specie] = new_table
    return new_ACE_table

def test():
    func = make_transition_function(0.7, 1.0, m=2)
    x = np.linspace(0.1, 3, 100)
    y = func(x)
    chik = compute_CHIK(x, 'SiO')
    func_zbl = make_ZBL_potential(14, 8) # SiO
    y_zbl = func_zbl(x)

    #ACE_path = '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/4/SiO2-4_24-20-16-12_pairpot.table'
    ACE_path = '/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-3-12REF_REP0.01CUT0.01_pairpot copy.table'
    ace_table = load_ACE_pairpot(ACE_path)
    print(ace_table['O_Si'])
    y_ACE = ace_table['O_Si'][:, 2]
    x_ACE = ace_table['O_Si'][:, 1]

    y_zbl_x_ACE = func_zbl(x_ACE)
    weight = func(x_ACE)
    mixed_func = y_ACE * weight + y_zbl_x_ACE * (1-weight)

    d_zbl = derivative(func_zbl, x_ACE, dx=1e-4)

    print(d_zbl)

    plt.plot(x,y)
    plt.savefig("transition_test.png")
    plt.clf()
    #plt.plot(x, chik, label='chik')
    plt.plot(x, y_zbl, label='zbl')
    plt.plot(x_ACE, y_ACE, label='ACE')
    plt.plot(x_ACE, mixed_func, linestyle = '--', color = 'red', label='mixed potential (ZBL + ACE)')
    plt.ylim((-20, 200))
    plt.xlim(right=3)
    plt.xlabel('r / Ang')
    plt.ylabel('V / eV')
    plt.legend()
    plt.savefig("chikzbltest.png")

def modify_ACE_pairpot_ZBL(path_potential):
    path_potential = Path(path_potential)
    savepath = path_potential.parent / (path_potential.stem+'_ZBL'+path_potential.suffix)
    print("saving ZBL modified pairpot table to " + str(savepath))
    ACE_table = load_ACE_pairpot(path_potential)
    ZBL_ACE_table = created_zbl_mixed_ACE(ACE_table)
    write_ACE_pairpot(savepath, ZBL_ACE_table)

if __name__ == '__main__':
    test()
    #exit()
    #path_potential = '/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-3-12REF_REP0.01CUT0.01_pairpot copy.table'
    julia_aces = {
        'ACE_12_0.5_0.0003' : '/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-12_s0.5_l0.0003_full',
        'ACE_18_N2750' : '/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-18_N2750',
        'ACE_REP1':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-12_s0.3_l0.0003_REP_full',
        'ACE_REP2':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-12_s0.3_l0.0003_REP2_full',
        'ACE_REP10':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-18_s0.3_l0.0003_REP10_full',
        'ACE_REP100':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-18_s0.3_l0.0003_REP100_full',
        'ACE_REP1000':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-18_s0.3_l0.0003_REP1000_full',
        'ACE_REP10000':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-18_s0.3_l0.0003_REP10000_full',
        'ACE_REP100000':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-18_s0.3_l0.0003_REP100000_full',
        'ACE_REP1000000':'/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-2-18_s0.3_l0.0003_REP1000000_full',
        'ACE_REP_0.01_REF': '/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-3-12REF_REP0.01',
        'ACE_REP_100_REF': '/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-3-12REF_REP100',
        'ACE_REP_0.01_REF_CUT0.01': '/mnt/scratch2/q13camb_scratch/adps2/summer_project/ace_fitting_julia/models/SiO2-3-12REF_REP0.01CUT0.01',

        'ACE_CHEAP_CHUCK_1': '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/1/SiO2-3-12',
        'ACE_CHEAP_CHUCK_2': '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/2/SiO2-3-18',
        'ACE_CHEAP_CHUCK_3': '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/3/SiO2-4-12',
        'ACE_CHEAP_CHUCK_4': '/mnt/scratch2/q13camb_scratch/adps2/ACE_NEW_CHUCK/4/SiO2-4_24-20-16-12',
        
        #'ACE_KAMIL_5': 'SiO2-2-8_s0.01_l0.01_full',
    }
    for p in julia_aces.values():
        modify_ACE_pairpot_ZBL(p+'_pairpot.table')
