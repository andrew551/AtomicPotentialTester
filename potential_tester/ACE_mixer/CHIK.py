import numpy as np
# from Carre, et al paper: New effective potential for Silica    (2008)
def compute_CHIK(r, specie):
    if specie == 'Si2':
        qq = 1.910418**2
        A = 3150.462646
        B = 2.851451
        C = 626.751953
        D = 3423200
    elif specie == 'SiO':
        qq = -1.910418**2 / 2
        A = 27029.419922
        B = 5.158606
        C = 148.099091
        D = 29
    elif specie == 'O2':
        qq = 1.910418**2 / 4
        A = 659.595398
        B = 2.590066
        C = 26.836679
        D = 113
    else:
        raise Exception("invalid specie")
    return qq * 14.3997 / r + A * np.exp(-B*r) - C / r**6 + D / r**24