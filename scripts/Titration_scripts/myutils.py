#!/usr/bin/python
import sys

AVOGADRO= 6.022140857E23

def get_nkaps_nfgs_boxsideA_at_desired_C_M(
        kap_C_M,
        fg_C_M,
        min_kap_number,
        min_fg_number,
        verbose=False):
    """
    Returns a tuple (n_kaps, n_fg_chains, box_side_A) that approximates kap_C_M and fg_C_M
    with at least min_kap_number and min_fg_number
    """
    global AVOGADRO
    if fg_C_M<1E-24:
        n_fg_chains= 0
    elif(kap_C_M>fg_C_M*10):
        n_fg_chains=min_fg_number
    else:
        n_fg_chains=max(100, min_fg_number)
    is_done=False
    while True:
      box_volume_L= n_fg_chains/AVOGADRO/fg_C_M
      box_volume_A3= box_volume_L*1E27
      box_side_A= box_volume_A3**(1/3.0)
      n_kaps_raw= box_volume_L*AVOGADRO*kap_C_M
      n_kaps = max(int(round(n_kaps_raw)), min_kap_number) if kap_C_M>0.0 else 0
      kap_C_M_actual= n_kaps/AVOGADRO/box_volume_L
      is_done= (kap_C_M==0.0 or \
                ( kap_C_M_actual/kap_C_M>0.5 and \
                  kap_C_M_actual/kap_C_M<2.0 ))
      if is_done:
         break
      else:
         n_fg_chains= n_fg_chains*2
         if verbose:
             print("Doubling n_fg_chains to {}".format(n_fg_chains))
    if verbose:
         print("n_kaps_raw: {}".format(n_kaps_raw))
         print "#n_fg_chains", n_fg_chains, "n_kaps", n_kaps
         print("#Kap concaentration %.3e [M]" % kap_C_M_actual)
         print("#FG concaentration %.3e [M]" % fg_C_M)
         print("#Molar ratio Kap:FG chains: %.3f" % (n_fg_chains/ (n_kaps+0.0)))
    return (n_kaps, n_fg_chains, box_side_A)


if __name__ == "__main__":
    sys.exit(0)
