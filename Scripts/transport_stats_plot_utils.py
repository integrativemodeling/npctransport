import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tempfile

def get_normalized_permeability(n_per_particle_per_sec, 
                                box_side_A = 2000,     
                                NE_thickness_A = 300,  
                                scale_factor = 1/100.0 
                                ):
    ''' c
    returns permeability in units of number of molecules per 
    second per NPC per uM concentration difference 
    
    :param n_per_particle_per_sec: number of transport events per particle per second simulation
    :param box_side_A: bounding box side length
    :param NE_thickness_A: nuclear enevelope thickness
    :param scale_factor: scaling factor for comparing to experimental data
    '''
    # Compute volume of one compartment 
    # (we want concentration difference between two sides of NPC per particle)
    vol_A3 = box_side_A**2 * (box_side_A-NE_thickness_A) / 2.0 
    vol_L = vol_A3 * 1e-27
    NA = 6.022e23
    C_particle_M = 1 / NA / vol_L
    C_particle_uM = C_particle_M * 1e6
    return n_per_particle_per_sec * scale_factor / C_particle_uM

def load_transport_stats(path, box_side_A = 2000, NE_thickness_A = 300,
                         scale_factor = 1/100.0):
    '''
    Load transport statistics from csv file and return a procecced dataframe.
    
    :param path: path to input csv file containing the following columns
        type: transport factor type (inert18 is an inert molecule with radius 18 A, kap18 is same for an NTR)
        n_per_particle_per_sec: number of transport events per particle per second simulation
        n_per_particle_per_sec_lbound: lower bound of n_per_particle_per_sec
        n_per_particle_per_sec_ubound: upper bound of n_per_particle_per_sec
    :param box_side_A: bounding box side length (for computing particle concentrations)
    :param NE_thickness_A: nuclear enevelope thickness (for computing particle concentrations)
    :param scale_factor: scaling factor for comparing to experimental data
    
    :return df: pandas dataframe with columns:
        is_kap: True if transport factor is a karyopherin
        radius: radius of transport factor in nm
        MW: molecular weight of transport factor in kDa
        n_per_sec_per_NPC_per_uM: permeability in units of number of molecules per
            second per NPC per uM concentration difference
        n_per_sec_per_NPC_per_uM_lbound: lower bound of permeability
        n_per_sec_per_NPC_per_uM_ubound: upper bound of permeability
    '''
    df = pd.read_csv(path)
    df['is_kap'] = df['type'].str.startswith('kap')
    df['radius'] = df['type'].str.extract(r'(\d+)').astype(int)
    df['MW'] = (df['radius'] / 20)**3 * 27
    df['n_per_sec_per_NPC_per_uM'] = \
        get_normalized_permeability(df['n_per_particle_per_sec'], 
                                    scale_factor=scale_factor)
    df['n_per_sec_per_NPC_per_uM_lbound'] = \
        get_normalized_permeability(df['n_per_particle_per_sec_lbound'], 
                                    scale_factor=scale_factor)
    df['n_per_sec_per_NPC_per_uM_ubound'] = \
        get_normalized_permeability(df['n_per_particle_per_sec_ubound'], 
                                    scale_factor=scale_factor)
    return df[['is_kap', 'radius', 'MW', 
               'n_per_sec_per_NPC_per_uM', 
               'n_per_sec_per_NPC_per_uM_lbound', 
               'n_per_sec_per_NPC_per_uM_ubound']]

def plot_permeability(ax, df, color, label, 
                      is_mw = True,
                      slope_location=None,  is_abbreviated_slope=False, is_slope_in_legend=False,
                      slope_width=0.5, slope_alpha=0.5,
                      verbose=False):
    '''
    plot NPC permeability for df with error bars using specified color and label
    If slope location is not None, add text with slope of fit in specified location
    
    :param ax: matplotlib axis to plot on
    :param df: pandas dataframe with columns:
        MW: molecular weight of transport factor in kDa (if is_mw is True)
        radius: molecular weight of transport factor in nm (if is_mw is False)
        n_per_sec_per_NPC_per_uM: permeability in units of number of molecules per
            second per NPC per uM concentration difference
        n_per_sec_per_NPC_per_uM_lbound: lower bound of permeability
        n_per_sec_per_NPC_per_uM_ubound: upper bound of permeability
    :label: label for plot
    :is_mw: if True, plot MW on x-axis, otherwise plot R
    :slope_location: location of text with slope of fit on normalized axes coordinates [0-1],[0-1]
                     if it is equal to 'last', plot text next to last data point
    :is_slope_in_legend: if True, add slope to legend labels, overriding slope_location
    :is_abbreviated_slope: if True, plot abbreviated slope text 
    :verbose: if True, print dataframe with permeability data 
    
    '''
    x = df['MW'] if is_mw else df['radius'] 
    y = df['n_per_sec_per_NPC_per_uM']
    y_err_low = df['n_per_sec_per_NPC_per_uM'] - df['n_per_sec_per_NPC_per_uM_lbound'] 
    y_err_high = df['n_per_sec_per_NPC_per_uM_ubound'] - df['n_per_sec_per_NPC_per_uM']
    y_is_positive = y > 0
    x = x[y_is_positive]
    y = y[y_is_positive]
    if verbose:
        print(label)
        display(pd.DataFrame({'R': np.round(20*(x/27)**0.3333), 'MW': x, 'p': y}))
        print("------------------")
    y_err_low = y_err_low[y_is_positive]
    y_err_high = y_err_high[y_is_positive]
    # add text with slope fitted to log(y) as a function of log(x), i.e. exponent of power law
    slope_dict = {}
    if (slope_location is not None or is_slope_in_legend) and len(x)>=2:
        log_x = np.log10(x)
        data_x = log_x if is_mw else x
        log_y = np.log10(y)
        data_y = log_y
        fit = np.polyfit(data_x, data_y, 1)
        fit_fn = np.poly1d(fit)
        ax.plot(x, 10**fit_fn(data_x),
                '-', color=color, 
                linewidth=slope_width, alpha=slope_alpha)
        if slope_location == 'last':
            transform = ax.transData
            ixmax = x.idxmax()
            slope_location = (x[ixmax]+0.1, y[ixmax])
        else:
            transform = ax.transAxes
        slope_dict[label] = r" $\propto " \
                            + ('MW' if is_mw else 'R') \
                            + "^{" \
                            + (f"{fit[0]:.1f}" if is_mw else f"{fit[0]:.2f}")  \
                            + r"}$"
        if is_slope_in_legend:
            label += slope_dict[label]
        else:
            ax.text(slope_location[0], slope_location[1],
                    (label if not is_abbreviated_slope else "")
                    + slope_dict[label],
                    transform=transform, 
                    color=color,
                    fontsize=15)
    # plot the main part:
    ax.errorbar(x, y, yerr=[y_err_low, y_err_high],
                fmt='o', color=color, label=label, alpha=0.75,
                markeredgecolor='none', 
                # errorbar style
                elinewidth=0.5, ecolor=color, capsize=1.5, capthick=0.5
    )


###### TESTS ######

def test_get_normalized_permeability():
    assert(np.isclose(get_normalized_permeability(3030.36, box_side_A = 2000, 
                                                  NE_thickness_A = 300,
                                                  scale_factor = 1/100.0),
                      62.04, atol=0.1))
    assert(np.isclose(get_normalized_permeability(3030.36, box_side_A = 2000, 
                                                  NE_thickness_A = 300,
                                                  scale_factor = 1.0),
                      6204, atol=1.0))

def test_load_transport_stats():
    CSV = b'''
From 10e-6,type,time_sec_mean,time_sec_sem,total_sim_time_sec,n_per_particle_per_sec,n_per_particle_per_sec_lbound,n_per_particle_per_sec_ubound,confidence_level
2,inert18,23.49993036472169,1.713908571425621,8.00002370569691,23.49993036472169,20.260645715762717,27.10983801762387,0.95
10,kap18,1407.620828920905,13.264693069640611,8.00002370569691,1407.620828920905,1381.7411038049772,1433.863494023603,0.95
'''
    with tempfile.NamedTemporaryFile() as fp:
        fp.write(CSV)
        fname = fp.name
        #fp.close()
        fp.seek(0)
        df = load_transport_stats(fname, box_side_A = 2000, NE_thickness_A = 300, scale_factor = 1/100.0)
    assert(len(df.index)==2)
    assert(len(df.columns)==6)
    assert(df['is_kap'].iloc[0] != df['is_kap'].iloc[1] == True)
    assert(df['radius'].iloc[0] == df['radius'].iloc[0] == 18)
    assert(np.isclose(df['MW'].iloc[0], 19.683, atol=0.01))
    assert(np.isclose(df['n_per_sec_per_NPC_per_uM'].iloc[0], 0.4812, atol=0.01))
    assert(np.isclose(df['n_per_sec_per_NPC_per_uM_lbound'].iloc[0], 0.4148, atol=0.01))
    assert(np.isclose(df['n_per_sec_per_NPC_per_uM_ubound'].iloc[0], 0.5550, atol=0.01))
        

def test_plot_utils():
    test_get_normalized_permeability()
    test_load_transport_stats()
    
test_plot_utils()