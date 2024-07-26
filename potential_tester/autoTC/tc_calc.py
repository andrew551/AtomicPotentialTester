import os

def calc_tc_hdf5_yaml(T_min,T_max,T_step,initial_string="phono3py_disp.yaml",output_file="tc.out",Q_MESH=[9,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py=None):
    
    if METHOD_thm :
        method="--thm"
    else:
        method= "--sigma=%f"%(sigma)+" --sigma-cutoff=%d"%(sigma_cutoff)


    settings= [ initial_string, 
                "--fc3 --fc2 --br --sym-fc",
                method,
                '--mesh="'+" ".join([str(s) for s in Q_MESH])+'"',
                "--isotope",
                "--tmin="+str(T_min),
                "--tmax="+str(T_max),
                "--tstep="+str(T_step),
                "--wigner"
        ]

    if path_phono3py is None:
        path_phono3py_folder=""
    else:
        path_phono3py_folder=path_phono3py+"/scripts/"
    
    path_phono3py_script=path_phono3py_folder+"phono3py "
    cmd=path_phono3py_script+" "+" ".join(settings)+" > "+output_file
    os.system('module load apps/python3/3.10.5/gcc-9.3.0 | python3 -m pip list | python3 -V')
    print(cmd, flush=True)
    os.system(cmd)

def calc_tc_hdf5_yaml_rta(T_min,T_max,T_step,initial_string="phono3py_disp.yaml",output_file="tc.out",Q_MESH=[9,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py=None):
    
    if METHOD_thm :
        method="--thm"
    else:
        method= "--sigma=%f"%(sigma)+" --sigma-cutoff=%d"%(sigma_cutoff)


    settings= [ initial_string, 
                "--fc3 --fc2 --br --sym-fc",
                method,
                '--mesh="'+" ".join([str(s) for s in Q_MESH])+'"',
                "--isotope",
                "--tmin="+str(T_min),
                "--tmax="+str(T_max),
                "--tstep="+str(T_step)
        ]

    if path_phono3py is None:
        path_phono3py_folder=""
    else:
        path_phono3py_folder=path_phono3py+"/scripts/"

    path_phono3py_script=path_phono3py_folder+"phono3py "
    cmd=path_phono3py_script+" "+" ".join(settings)+" > "+output_file
    print(cmd)

    os.system(cmd)

def calc_tc_params(T_min,T_max,T_step,initial_string="phono3py_params.yaml",output_file="tc.out",Q_MESH=[9,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py=None):
    
    if METHOD_thm :
        method="--thm"
    else:
        method= "--sigma=%f"%(sigma)+" --sigma-cutoff=%d"%(sigma_cutoff)


    settings= [ initial_string, 
                "--br --sym-fc",
                method,
                '--mesh="'+" ".join([str(s) for s in Q_MESH])+'"',
                "--isotope",
                "--tmin="+str(T_min),
                "--tmax="+str(T_max),
                "--tstep="+str(T_step),
                "--wigner",
                "--nonac"
        ]

    if path_phono3py is None:
        path_phono3py_folder=""
    else:
        path_phono3py_folder=path_phono3py+"/scripts/"

    path_phono3py_script=path_phono3py_folder+"phono3py-load "
    cmd=path_phono3py_script+" "+" ".join(settings)+" > "+output_file
    print(cmd)

    os.system(cmd)


def calc_tc_params_nac(T_min,T_max,T_step,initial_string="phono3py_params.yaml",output_file="tc.out",Q_MESH=[11,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py=None):

    if METHOD_thm :
        method="--thm"
    else:
        method= "--sigma=%f"%(sigma)+" --sigma-cutoff=%d"%(sigma_cutoff)


    settings= [ initial_string,
                "--br --sym-fc",
                method,
                '--mesh="'+" ".join([str(s) for s in Q_MESH])+'"',
                "--isotope",
                "--tmin="+str(T_min),
                "--tmax="+str(T_max),
                "--tstep="+str(T_step),
                "--wigner"
        ]

    if path_phono3py is None:
        path_phono3py_folder=""
    else:
        path_phono3py_folder=path_phono3py+"/scripts/"

    path_phono3py_script=path_phono3py_folder+"phono3py-load "
    cmd=path_phono3py_script+" "+" ".join(settings)+" > "+output_file
    print(cmd)

    os.system(cmd)


def calc_tc_params_nac_rta(T_min,T_max,T_step,initial_string="phono3py_params.yaml",output_file="tc.out",Q_MESH=[11,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py=None):

    if METHOD_thm :
        method="--thm"
    else:
        method= "--sigma=%f"%(sigma)+" --sigma-cutoff=%d"%(sigma_cutoff)


    settings= [ initial_string,
                "--br --sym-fc",
                method,
                '--mesh="'+" ".join([str(s) for s in Q_MESH])+'"',
                "--isotope",
                "--tmin="+str(T_min),
                "--tmax="+str(T_max),
                "--tstep="+str(T_step)
        ]

    if path_phono3py is None:
        path_phono3py_folder=""
    else:
        path_phono3py_folder=path_phono3py+"/scripts/"

    path_phono3py_script=path_phono3py_folder+"phono3py-load "
    cmd=path_phono3py_script+" "+" ".join(settings)+" > "+output_file
    print(cmd)

    os.system(cmd)

def calc_tc_disp_hdf5_nac_rta(T_min,T_max,T_step,initial_string="phono3py.yaml",output_file="tc.out",Q_MESH=[9,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py=None,output_hdf5=None):

    if METHOD_thm :
        method="--thm"
    else:
        method= "--sigma=%f"%(sigma)+" --sigma-cutoff=%d"%(sigma_cutoff)


    settings= [ initial_string,
                "--br --fc2 --fc3",
                method,
                '--mesh="'+" ".join([str(s) for s in Q_MESH])+'"',
                "--isotope",
                "--tmin="+str(T_min),
                "--tmax="+str(T_max),
                "--tstep="+str(T_step),
                "--nac",
                "-o "+output_hdf5 if output_hdf5 is not None else ""
        ]

    if path_phono3py is None:
        path_phono3py_folder=""
    else:
        path_phono3py_folder=path_phono3py+"/scripts/"

    path_phono3py_script=path_phono3py_folder+"phono3py "
    cmd=path_phono3py_script+" "+" ".join(settings)+" > "+output_file
    print(cmd)

    os.system(cmd)

def calc_tc_disp_hdf5_rta(T_min,T_max,T_step,initial_string="phono3py.yaml",output_file="tc.out",Q_MESH=[9,9,9],METHOD_thm=True,sigma=0.012,sigma_cutoff=3,path_phono3py=None,output_hdf5=None):

    if METHOD_thm :
        method="--thm"
    else:
        method= "--sigma=%f"%(sigma)+" --sigma-cutoff=%d"%(sigma_cutoff)


    settings= [ initial_string,
                "--br --fc2 --fc3",
                method,
                '--mesh="'+" ".join([str(s) for s in Q_MESH])+'"',
                "--isotope",
                "--tmin="+str(T_min),
                "--tmax="+str(T_max),
                "--tstep="+str(T_step),
                "-o "+output_hdf5 if output_hdf5 is not None else ""
        ]

    if path_phono3py is None:
        path_phono3py_folder=""
    else:
        path_phono3py_folder=path_phono3py+"/scripts/"

    path_phono3py_script=path_phono3py_folder+"phono3py "
    cmd=path_phono3py_script+" "+" ".join(settings)+" > "+output_file
    print(cmd)

    os.system(cmd)