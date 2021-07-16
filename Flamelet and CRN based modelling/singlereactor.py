import cantera as ct
import math
import sys
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14, 'legend.fontsize':14, 'lines.linewidth':3})
from cycler import cycler
#plt.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
import numpy
from numpy import dot, array, add, matrix, linalg, subtract, ndarray, float64,zeros, allclose, matmul, multiply, transpose, append, amax, amin
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.linalg import spsolve, gmres, lgmres, splu, spilu, inv

gas=ct.Solution('gri30.cti')
# species taken from gri-3.0
species = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'C', 'CH', \
           'CH2', 'CH2(S)', 'CH3', 'CH4', 'CO', 'CO2', 'HCO', 'CH2O', \
           'CH2OH', 'CH3O', 'CH3OH', 'C2H', 'C2H2', 'C2H3', 'C2H4', 'C2H5', \
           'C2H6', 'HCCO', 'CH2CO', 'HCCOH', 'N', 'NH', 'NH2', 'NH3', 'NNH', \
           'NO', 'NO2', 'N2O', 'HNO', 'CN', 'HCN', 'H2CN', 'HCNN', 'HCNO', 'HOCN', \
           'HNCO', 'NCO', 'N2', 'AR', 'C3H7', 'C3H8', 'CH2CHO', 'CH3CHO']

T=1500
fig , ax = plt.subplots()
ax.set_xlabel('Time(s)',fontweight='bold')
#ax1=ax.twinx()
col=['tab:cyan','tab:red','tab:blue','tab:green']
ener='on'
a=[1]
b=[2,3]
c=a+b
print max([1,2])
y= gas.Y

a=numpy.array([1,1,1])
b=numpy.array([2,2,2])
a=-1*a
print a
c=numpy.add(a,b)
print c
R = ct.Reaction.listFromFile('gri30.cti')
mat=[[0]*len(y)]*len(R)
for sp in range(len(y)):
    y0=y
    dy=0.01*y0[sp]
    y0[sp]=y0[sp]+dy
    r0=gas.net_rates_of_progress
    gas.Y=y0
    r=gas.net_rates_of_progress
    for rat in range(len(r)):
        mat[rat][sp]=(r[rat]-r0[rat])/(dy+1e-12)
mat_mu=[[0]*len(y)]*len(R)
print r
for r in range(len(R)):
    react=R[r].reactants
    for s in react.keys():
        ind= gas.species_index(s)
        mat_mu[r][ind]=-react[s]
    prod=R[r].products
    for s in prod.keys():
        ind = gas.species_index(s)
        mat_mu[r][ind] = prod[s]
print react
print prod
Jw=numpy.zeros((len(y),len(y)))
for ind in range(len(R)):
    S=numpy.matrix(mat_mu[ind])
    print "S=",S
    dR=numpy.matrix([mat[ind]])
    print "dr=",dR
    S=numpy.transpose(S)
    Jw_temp=numpy.matmul(S,dR)
    print numpy.size(Jw_temp)
    Jw=numpy.add(Jw,Jw_temp)

print Jw
print Jw[(1,1)]


for temp in range(3):
    vol=0.00039519
    num=1
    PSR=[]
    m1=2.56e-5
    m2=0.00057425
    m3=m1+m2

    for i in range(num):
        gas.TPY= T, ct.one_atm, {'O2':0.23,'N2':0.77}
        psr=ct.IdealGasReactor(gas, energy=ener)
        if i==0:
            psr.chemistry_enabled=True
        psr.volume=vol/num
        PSR.append(psr)
        if i!=0:
            mfc_con=ct.MassFlowController(PSR[i-1],PSR[i],mdot=m3)
    rates1=psr.thermo.net_production_rates
    conc=psr.thermo.concentrations
    print PSR[0].mass
    R=ct.Reaction.listFromFile('gri30.cti')
    print "dir=",PSR[0].thermo.reverse_rates_of_progress
    print "orders=",R[5].reactants

    #print conc
    moles=(psr.thermo.P*psr.volume)/(psr.thermo.T*83140)
    tscale=[]
    for r in range(len(rates1)):
        tscale.append(conc[r]/((rates1[r])+1e-12))
        #tscale.append(abs((conc[r]**2)/((rates1[r])+1e-12)))
    print "tscale_max=",max(tscale)
    print "tscale_min=",min(tscale)

    gas.TPY=300, ct.one_atm, {'CH4':1}
    print gas.T
    res1=ct.Reservoir(gas)
    gas.TPX=673, ct.one_atm, {'O2':0.21,'N2':0.79}
    res2=ct.Reservoir(gas)
    res3=ct.Reservoir(gas)
    print res2.thermo.T
    mfc1=ct.MassFlowController(res1,PSR[0],mdot=m1)
    mfc2=ct.MassFlowController(res2,PSR[0],mdot=m2)
    mfc3=ct.MassFlowController(PSR[num-1],res3,mdot=m3)
    kv = PSR[num - 1].volume
    #kv=(PSR[num-1].mass/m3)*((PSR[num-1].volume)**(1/3))
    print "kv=",kv
    print mfc1.mdot(0)
    kv1=m1
    kv2=m2
    kv3=m3
    valve0 = ct.Valve(PSR[0], res1)
    valve0.set_valve_coeff(kv1)
    valve01 = ct.Valve(PSR[0], res2)
    valve01.set_valve_coeff(kv2)
    valve02 = ct.Valve(res3,PSR[num-1])
    valve02.set_valve_coeff(kv3)
    valve1=ct.Valve(PSR[num-1],res3)
    valve1.set_valve_coeff(kv3)#lambda dP: (1e-5 * dP) ** 2)
    print "valve dir=", dir(valve0)
    net=ct.ReactorNet(PSR)
    xin= PSR[0].thermo.Y
    print "xin=",xin
    net.advance_to_steady_state()
    print "xin fin=",xin
    gas.TPX = 300.0, ct.one_atm, 'H:1.0'
    """igniter = ct.Reservoir(gas)
    # The igniter will use a Gaussian time-dependent mass flow rate.
    fwhm = 0.2
    amplitude = 0.00001
    t0 = 0.0001
    igniter_mdot = lambda t: amplitude * math.exp(-(t - t0) ** 2 * 4 * math.log(2) / fwhm ** 2)
    m_ignite = ct.MassFlowController(igniter, PSR[0], mdot=igniter_mdot)"""

    gas.TPY=PSR[0].thermo.T,PSR[0].thermo.P,PSR[0].thermo.Y
    y=PSR[0].thermo.Y
    R = ct.Reaction.listFromFile('gri30.cti')
    # Initializing matrix of dimension NrxNs
    mat = [[0] * len(y)] * len(R)
    # Assembling matrix of dR/dY, where each row is for a particular reaction
    r0 = gas.net_rates_of_progress
    y0 = y
    for sp in range(len(y)):


        dy = 0.01
        y0[sp] = y0[sp] + dy
        print y0
        gas.TPY =PSR[0].thermo.T,PSR[0].thermo.P, y0
        r = gas.net_rates_of_progress
        y0[sp] = y0[sp] - dy
        for rat in range(len(r)):
            mat[rat][sp] = (r[rat] - r0[rat]) / (dy + 1e-12)
    mat_mu = [[0] * len(y)] * len(R)
    for r in range(len(R)):
        # dictionary of reactant stoichiometric coeff
        react = R[r].reactants
        for s in react.keys():
            ind = gas.species_index(s)
            mat_mu[r][ind] = -react[s] * gas.molecular_weights[ind]
        # dictionary of product sotichiometric coeff
        prod = R[r].products
        for s in prod.keys():
            ind = gas.species_index(s)
            mat_mu[r][ind] = prod[s] * gas.molecular_weights[ind]
    Jw = zeros((len(y), len(y)))
    for ind in range(len(R)):
        S = matrix(mat_mu[ind])
        dR = matrix(mat[ind])
        S = transpose(S)
        Jw_temp = matmul(S, dR)
        Jw = add(Jw, Jw_temp)
    print "Jw=", Jw
    #globalsolution
    Jac = []
    row_Jac = []
    col_Jac = []
    Cw = [0] * (len(species) * len(PSR))
    recalc = 1
    rhsvect=[m3]
    RatesofProduction = [0] * (len(species) * len(PSR))
    f_vect = [0] * (len(species) * len(PSR))
    # Cw = [0] * (len(species) * len(PSR))
    SpecMassFrac = [0] * (len(species) * len(PSR))
    # Obtaining rates of production
    for r_ind in range(len(PSR)):
        r = PSR[r_ind]
        xfin = r.thermo.Y
        rates = r.thermo.net_production_rates
        for ind in range(len(rates)):
            RatesofProduction[r_ind * len(species) + ind] = rates[ind] * r.thermo.molecular_weights[ind] / r.density
            f_vect[r_ind * len(species) + ind] = rhsvect[r_ind] * xfin[ind] / r.mass
            SpecMassFrac[r_ind * len(species) + ind] = xfin[ind]
            for pt in range(len(PSR)):
                Cw[r_ind * len(species) + ind] += -1 * m3 * PSR[pt].thermo.Y[
                    ind] / r.mass
    RHS_Jac = add(array(Cw), array(f_vect))
    RHS_Jac = add(RHS_Jac, array(RatesofProduction))
    RHS_Jac = -1 * RHS_Jac
    print "Jacobian Calculation"
    for p in range(len(PSR)):
        if recalc == 0:
            q = p
            coeff = 1
            for spec0 in range(len(species)):
                for spec in range(len(species)):
                    Js = 0
                    Jw_add = 0
                    if spec == spec0:
                        Js = -1 * coeff * m3 / PSR[p].mass

                    # Column in source term jacobian has constant denominator species
                    Jw_add = Jw[spec, spec0] / PSR[p].density
                    J = Js + Jw_add
                    if J != 0:
                        # Row in jacobian has constant denominator species
                        v1 = q * len(species) + spec0
                        v2 = p * len(species) + spec
                        r_index = row_Jac.index(v1)
                        c_index = col_Jac.index(v2)
                        ind_start = max([r_index, c_index])
                        for row_ind in range(len(row_Jac)):
                            ind_update = row_ind + ind_start
                            try:
                                if row_Jac[ind_update] == v1 and col_Jac[ind_update] == v2:
                                    Jac[ind_update] = J
                                    break
                            except:
                                Jac.append(J)
                                row_Jac.append(v1)
                                col_Jac.append(v2)


        else:
            for q in range(len(PSR)):
                coeff = 1
                if coeff != 0:
                    for spec0 in range(len(species)):
                        for spec in range(len(species)):
                            Js = 0
                            Jw_add = 0
                            if spec == spec0:
                                Js = -1 * coeff * m3 / PSR[p].mass
                            if p == q:
                                # Column in source term jacobian has constant denominator species
                                Jw_add = Jw[spec, spec0] / PSR[p].density
                            J = Js + Jw_add
                            if J != 0:
                                Jac.append(J)
                                # Row in jacobian has constant denominator species
                                row = p * len(species) + spec
                                col = q * len(species) + spec0
                                row_Jac.append(row)
                                col_Jac.append(col)

    PSR_Jac = csc_matrix((Jac, (row_Jac, col_Jac)))
    """print "Preconditioning"
    Precond=splu(PSR_Jac)
    SpecMassFrac_update = Precond.solve(RHS_Jac)
    print "Inverting L and U"
    u_1=inv(Precond.U)
    l_1=inv(Precond.L)
    print "Computing M"
    M=u_1*l_1"""
    print "Starting LGMRES"
    SpecMassFrac_update, exitcode = lgmres(A=PSR_Jac, b=RHS_Jac, maxiter=40)
    print "Exitcode=", exitcode
    print SpecMassFrac_update
    e_list = []
    for e in range(len(SpecMassFrac_update)):
        e_list.append(abs(SpecMassFrac_update[e] / (SpecMassFrac[e] + 1e-12)))
    e_newt_max = amax(e_list)
    error = e_newt_max
    SpecMassFrac_update2 = add(SpecMassFrac_update, SpecMassFrac)
    print SpecMassFrac_update2
    print "error=", error
    for p in PSR:
        begin = PSR.index(p) * len(species)
        y0 = []
        for i in range(len(species)):
            y0.append(SpecMassFrac_update2[begin + i])
        g = p.thermo
        Temp = g.T
        Pr = g.P
        g.TPY = Temp, Pr, y0
        p.insert(g)

    print "PSR op"
    print PSR[0].thermo.Y
    time=[]#[0]
    t=0
    list_spec=psr.thermo.species_names
    index=list_spec.index('CO2')
    index2=list_spec.index('NO')
    index3 = list_spec.index('O')
    Temp=[]#[PSR[num-1].thermo.T]
    conc=[]#[PSR[num-1].thermo.Y[index]]
    conc2=[]#[PSR[num-1].thermo.Y[index2]]
    conc3 = []#[PSR[num - 1].thermo.Y[index3]]
    rates2 = psr.thermo.net_rates_of_progress
    reac_num=177
    #reac_num-=1
    nox_rate=[]#[rates2[reac_num]]
    tres=[PSR[num-1].mass/m3]
    print tres
    dt=0.0001#m3/psr.mass
    tstart=0
    error=0
    e=[]
    max_it=20000
    error_log=[]
    #print PSR[0].thermo.reactions()
    reactionoff=[177,183,184,203,204,205,238,239,241]
    for i in range(4):
        reactionoff.append(177+i)
    for i in range(5):
        reactionoff.append(185+i)
    for i in range(5):
        reactionoff.append(211+i)
    for i in range(15):
        reactionoff.append(243+i)
    reactionoff.append(194)
    reactionoff.append(197)
    reactionoff.append(198)
    reactionoff.append(207)
    reactionoff.append(221)
    reactionoff.append(223)
    reactionoff.append(225)
    reactionoff.append(227)
    reactionoff.append(228)
    reactionoff.append(273)
    reactionoff.append(280)
    reactionoff.append(282)




    for i in range(max_it):
        inp=input("Enter something")
        if inp=='stop':
            sys.exit()
        if T > 100:
            for r in reactionoff:
                PSR[0].kinetics.set_multiplier(0.0, r)
     
        net.set_initial_time(tstart)
        xin=PSR[num-1].thermo.Y
        net.advance(tstart+dt)

        #net.advance_to_steady_state()
        xfin=PSR[num-1].thermo.Y
        if i==max_it-1:
            for y in range(len(xfin)):
                e.append(abs(xfin[y]-xin[y])/(xin[y]+1e-15))
            error=max(e)
            error_log.append(error)
            max_ind=e.index(error)
        tstart+=dt
        #net.advance_to_steady_state()
        time.append(tstart)
        t=t+dt
        rates_nox = -psr.thermo.net_rates_of_progress[reac_num]
        nox_rate.append(rates_nox)
        Temp.append(PSR[num-1].thermo.T)
        conc.append(PSR[num-1].thermo.Y[index])
        conc2.append(PSR[num-1].thermo.Y[index2])
        conc3.append(PSR[num - 1].thermo.Y[index3])
        tres.append(PSR[num-1].mass/m3)
    print 'error=',error
    #h=ax.plot(time, conc, label=list_spec[index]+' %iK'%T, color=col[temp])
    if temp>100:
        h = ax.plot(time, Temp, label='T0=%iK, NO off' % T, color=col[temp])
        #h = ax.plot(time, conc, label=list_spec[index] + ' @T0=%iK, NO off' % T, color=col[temp])
    else:
        h = ax.plot(time, Temp, label='T0=%iK' % T, color=col[temp])
        #h = ax.plot(time, conc, label=list_spec[index] + ' @T0=%iK' % T, color=col[temp])
    cid=h[0].get_color()
    ax.set_ylabel('Temperature', fontweight='bold')
    #ax.set_ylabel('Pressure(Pa)', fontweight='bold')
    #ax.set_ylabel('Mass fraction of '+list_spec[index],fontweight='bold')
    #ax1.plot(time, nox_rate, label='Thermal NO reaction rate'+'(%iK)' % T, color=cid, linestyle='--')
    """if temp>0:
        #h = ax.plot(time, Temp, label='T0=%iK, NO off' % T, color=col[temp])
        ax1.plot(time, conc2, label=list_spec[index2] + ' @T0=%iK, NO off' % T, color=col[temp], linestyle='--')
    else:
        #h = ax.plot(time, Temp, label='T0=%iK' % T, color=col[temp])
        ax1.plot(time, conc2, label=list_spec[index2] + ' @T0=%iK' % T, color=col[temp], linestyle='--')
    #ax1.plot(time, conc2, label=list_spec[index2]+' @T0=%iK, NO off'%T, color=cid, linestyle='--')
    ax1.set_ylabel('Mass fraction of '+list_spec[index2],fontweight='bold')"""
    #ax1.set_ylabel('Net rate of progress(kmol/(m^3-s))', fontweight='bold')
    #ax.plot(time, conc3, label=list_spec[index3]+' %iK' % T, color=cid, linestyle='-.')
    #ax.set_ylabel(list_spec[index3])
    T=T+600
    #print "error co2=",e[list_spec.index('CO2')]
    print "error=",error
    #print "Species=",list_spec[max_ind]
    #print "tres=",tres
    for i in PSR:
        print i.thermo.T
    
h0,lbl0= ax.get_legend_handles_labels()
#h1,lbl1= ax1.get_legend_handles_labels()
h=h0#+h1
lbl=lbl0#+lbl1
plt.legend(h,lbl,loc=4)
#plt.savefig('singlereactor_no_prod.png')
plt.show()   
