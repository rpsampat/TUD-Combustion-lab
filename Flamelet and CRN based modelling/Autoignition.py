import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
gas = ct.Solution('POLIMI_C1C3.cti')
gas.transport_model = 'Multi'
print dir(ct.Solution)
#gas = ct.Solution('POLIMI_C1C3.cti')
phi = 0.0
phi_list = []
Temp_auto = []
Time_auto = []
dt = 1e-5
flowrate=0.001
T=973
fig1, ax1 = plt.subplots()
# plt.yscale('log')
ax1.set_ylabel('T$_{mixture}$(K)')
ax1.set_xlabel('Fresh reactant concentration')
#ax1.xaxis.set_major_locator(MultipleLocator(0.1))
#ax1.xaxis.set_minor_locator(MultipleLocator(0.02))
ax1.grid(True, ls='--')
#equivalence ratio variation of mixture
for i in range(5):
    phi = 0.5#phi + 0.4
    #T = T + 100  # 1000.0
    a = 2.0 / phi
    o2_perc = 0.21-(i*0.04)
    b = a*((1-o2_perc)/o2_perc)-1
    phi_list.append(o2_perc*100)
    MWair=28.85
    MWfuel=16.0
    m_air=4.76*a*MWair
    m_fuel=1.0*MWfuel
    m_tot=m_air+m_fuel
    Y_fuel=m_fuel/m_tot
    Y_air=m_air/m_tot

    Temp_auto.append(10.0)
    Time_auto.append(10.0)
    #temperature variation of reactants
    for j in range(1):
        print 'subiteration %i'%(j)
        #T = T + 100.0

        gas.TPX = T, ct.one_atm, {'CH4': 1.0, 'O2': a, 'N2': b}
        Cp=gas.cp
        print "Cp=",Cp
        res=ct.Reservoir(gas)
        res2=ct.Reservoir(gas)
        reactor = ct.IdealGasConstPressureReactor(gas, energy='on')
        #reactor = ct.IdealGasReactor(gas, energy='on')
        print reactor.thermo.thermal_conductivity
        # mfc1 = ct.MassFlowController(res, reactor)
        # mfc1.set_mass_flow_rate(flowrate)
        # mfc2 = ct.MassFlowController(reactor, res2)
        # mfc2.set_mass_flow_rate(flowrate)
        net = ct.ReactorNet([reactor])
        #res = ct.Reservoir(gas)
        #res2= ct.Reservoir(gas)
        #v1=ct.Valve(reactor, res)
        #v1.set_valve_coeff(reactor.volume)
        #m1 = ct.MassFlowController(res,reactor,mdot=1)
        #m2 = ct.MassFlowController(reactor, res2, mdot=1e-2)
        #pc = ct. PressureController(reactor, res, master=m1, K=reactor.volume)
        t = 0.0
        # try:
        #     net.advance_to_steady_state()
        # except:
        #     pass
        T_temp=[]
        time_temp=[]

        for iter in range(1000000):
            t = t+dt
            net.advance(t)
            T_temp.append(reactor.thermo.T)
            time_temp.append(t)
        T_deriv=list(np.diff(T_temp))
        ax1.plot(time_temp, T_temp)
        print np.max(T_temp)
        ind=T_deriv.index(np.max(T_deriv))
        print ind
        t_print=time_temp[ind-1]
        T_final = reactor.thermo.T
        print T_final
        Temp_auto[len(Temp_auto) - 1] = T_final#t_print#
        Time_auto[len(Time_auto) - 1] = t_print#
        if T_final - T > 100.0:
            #Temp_auto[len(Temp_auto)-1]=T
            print T_final
            break




fig2 , ax2 = plt.subplots()
ax2.plot(phi_list, Temp_auto,'*')

fig3 , ax3 = plt.subplots()
ax3.plot(phi_list, Time_auto,'*')
plt.show()
