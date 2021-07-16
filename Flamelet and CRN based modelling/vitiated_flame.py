import cantera as ct
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import numpy as np


def vitiated_product(T, phi, Tset):
    gas = ct.Solution('gri30.cti')
    species = gas.species_names
    a = 2.0/phi
    gas.TPX = T, ct.one_atm, {'CH4': 1.0, 'O2': a, 'N2': a*3.76}#{'CO2':0.9,'OH':0.1}#{'CO2':0.5,'O2':0.5}#
    #reactor = ct.IdealGasReactor(gas, energy='on')
    reactor = ct.IdealGasConstPressureReactor(gas, energy='on')
    net = ct.ReactorNet([reactor])
    net.advance_to_steady_state()
    Tout = gas.T
    Xco2 = 0#reactor.thermo.X[species.index('CO2')]
    Xn2 = reactor.thermo.X[species.index('N2')]
    Xh2o = 0#reactor.thermo.X[species.index('H2O')]
    Xoh = reactor.thermo.X[species.index('OH')]
    Xh = reactor.thermo.X[species.index('H')]
    gas.TPX = Tset, ct.one_atm, gas.X#{'CO2': Xco2, 'H2O': Xh2o, 'N2': Xn2, 'OH':Xoh, 'H':Xh}
    return gas


def fresh_reactant(T, phi):
    a = 2.0 / phi
    gas = ct.Solution('gri30.cti')
    gas.TPX = T, ct.one_atm, {'CH4': 1.0, 'O2': a, 'N2': a * 3.76}
    return gas


def gas_mixture(gas1, gas2, option, T_mix, compo, plot_enabled):
    dt=1e-5
    epsilon = 1e-5
    species = gas1.species_names
    h1, p1 = gas1.HP
    h2, p2 = gas2.HP
    time_autoign = []
    time_char = []
    T_rise = []
    T_initial_list=[]
    if plot_enabled == 'y':
        fig1, ax = plt.subplots()
        fig_ch4, ax_ch4 = plt.subplots()
    for Y in compo:
        ch4_conc = []
        h_conc = []
        T_temp = []
        time_temp = [0]
        h=((1-Y)*h1+Y*h2)
        compo_gas2=gas2.Y
        compo_mix=gas1.Y
        for i in range(len(species)):
            compo_mix[i]=(1-Y)*compo_mix[i]+Y*compo_gas2[i]
        gas = ct.Solution('gri30.cti')
        if option=='A':
            gas.TPY = T_mix, ct.one_atm, compo_mix
        else:
            gas.HPY = h, ct.one_atm, compo_mix
        T_initial=gas.T
        print 'Tinitial=',T_initial
        # Defining reactor
        reactor = ct.IdealGasConstPressureReactor(gas, energy='on')
        #reactor = ct.IdealGasReactor(gas)
        net = ct.ReactorNet([reactor])
        # Solving reactor by time stepping
        t = 0.0
        T_temp.append(reactor.thermo.T)
        ch4_conc.append(reactor.thermo.Y[species.index('CH4')])
        h_conc.append(reactor.thermo.Y[species.index('H')])
        for iter in range(100000):
            T_temp.append(reactor.thermo.T)
            time_temp.append(t)
            ch4_conc.append(reactor.thermo.Y[species.index('CH4')])
            h_conc.append(reactor.thermo.Y[species.index('H')])
            net.advance(t)
            t = t + dt
            # Post processing quantities
        T_final = reactor.thermo.T
        #print 'Tfinal=', T_final
        T_deriv = list(np.diff(T_temp))
        H_deriv = list((h_conc))
        delta_T=T_final-T_initial
        print 'Delta T=', delta_T
        # Determining autoignition delay time and chemical timescale
        max_T_deriv=np.max(T_deriv)
        #ind = T_deriv.index(max_T_deriv)
        max_H_deriv = np.max(H_deriv)
        ind = H_deriv.index(max_H_deriv)
        T_deriv_sub = T_deriv[ind:]
        for ind_steady in range(len(T_deriv_sub)):
            if T_deriv_sub[ind_steady]<epsilon:
                break
        # Characteristic timescale defined as time from start of reaction upto steady state of temperature
        try:
            t_char=time_temp[ind+ind_steady]*1e3
        except:
            t_char = time_temp[iter] * 1e3
        # Autoignition delay time defined by time at which maximum derivative of temperature occurs
        if max_T_deriv<=1e-3 or delta_T<2:
            t_auto=1
        else:
            t_auto = time_temp[ind - 1]

        t_auto = time_temp[ind - 1]

        # Read out quantities
        if plot_enabled=='y':
            ax.plot(time_temp, T_temp)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Temperature (K)')

            #ax_ch4.plot( time_temp, ch4_conc)
            ax_ch4.plot(ch4_conc)
            ax_ch4.set_xlabel('Time (s)')
            ax_ch4.set_ylabel('Y$_{CH4}$')

        time_autoign.append(t_auto*1e6)
        time_char.append(t_char)
        T_rise.append(delta_T)
        T_initial_list.append(T_initial)

    return time_autoign, time_char, T_rise,T_initial_list


if __name__=='__main__':
    # Vitiated product
    phi_prod = 0.5
    T_prod = 1100  # K
    gas1 = vitiated_product(T_prod, phi_prod, T_prod)
    print gas1.T

    # Fresh reactant
    phi = 0.9
    T= 1300  # K
    gas2 = fresh_reactant(T, phi)

    # Figures
    """fig2, ax2 = plt.subplots()
    #ax2.set_yscale('log')
    ax2.set_ylabel('Autoignition-delay time ($\mu$s)')
    ax2.set_xlabel('Fresh reactant concentration')
    ax2.xaxis.set_major_locator(MultipleLocator(0.1))
    ax2.xaxis.set_minor_locator(MultipleLocator(0.02))
    #ax2.set_yticks([50, 100, 500, 1000, 2000])
    #ax2.set_yticklabels([50, 100, 500, 1000, 2000])

    ax2.grid(True, ls='--')

    fig3, ax3 = plt.subplots()
    # plt.yscale('log')
    ax3.set_ylabel('Characteristic Chemical Time scale(ms)')
    ax3.set_xlabel('Fresh reactant concentration')
    ax3.xaxis.set_major_locator(MultipleLocator(0.1))
    ax3.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax3.grid(True, ls='--')

    fig4, ax4 = plt.subplots()
    # plt.yscale('log')
    ax4.set_ylabel('$\Delta$T(K)')
    ax4.set_xlabel('Fresh reactant concentration')
    ax4.xaxis.set_major_locator(MultipleLocator(0.1))
    ax4.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax4.grid(True, ls='--')

    fig5, ax5 = plt.subplots()
    # plt.yscale('log')
    ax5.set_ylabel('T$_{mixture}$(K)')
    ax5.set_xlabel('Fresh reactant concentration')
    ax5.xaxis.set_major_locator(MultipleLocator(0.1))
    ax5.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax5.grid(True, ls='--')"""

    # Mixture
    compo = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9, 0.98]#[0.01, 0.07, 0.08, 0.1,
    option = 'B' # 'A': T_mix used as starting temperature, 'B': enthalpy balance based on mixing of 2 gases
    T_mix = [1600]#[1700, 1900, 2000] # K
    col=['k','b','g','r']
    Tstart = 1200
    for i in range(2):

        t_mix = Tstart + i*50
        gas1 = vitiated_product(T_prod, phi_prod, t_mix)
        time_autoign, time_char, T_rise, T_initial = gas_mixture(gas1, gas2, option, t_mix, compo, plot_enabled='n')
        try:
            time_autoign_arr = np.vstack((time_autoign_arr, time_autoign))
            time_char_arr = np.vstack((time_char_arr, time_char))
            T_rise_arr = np.vstack((T_rise_arr, T_rise))
        except:
            time_autoign_arr = np.array(time_autoign)
            time_char_arr = np.array(time_char)
            T_rise_arr = np.array(T_rise)
        """ax2.plot(compo, time_autoign,'*', markerfacecolor=col[T_mix.index(t_mix)], mec='none', ms = 8)
        ax3.plot(compo, time_char,'o', markerfacecolor=col[T_mix.index(t_mix)], mec='none', ms = 8)
        ax4.plot(compo, T_rise,'*', markerfacecolor=col[T_mix.index(t_mix)], mec='none', ms =8)
        ax5.plot(compo, T_initial, '*', markerfacecolor=col[T_mix.index(t_mix)], mec='none', ms=8)"""
    '''T_labels=[str(i)+'K' for i in T_mix]
    ax2.legend(T_labels)
    ax3.legend(T_labels)
    ax4.legend(T_labels)'''
    extent = min(compo), max(compo), Tstart, t_mix
    fig6 = plt.figure()
    imgplot6 = plt.imshow(time_autoign_arr, cmap='seismic', origin="lower")
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Auto ignition delay")
    plt.xlabel('Fresh reactant fraction')
    plt.ylabel('Tstart')

    fig7 = plt.figure()
    imgplot7 = plt.imshow(time_char_arr, cmap='seismic')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Characteristic time")
    fig8 = plt.figure()
    imgplot8 = plt.imshow(T_rise_arr, cmap='seismic')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Temperature rise")

    plt.show()



