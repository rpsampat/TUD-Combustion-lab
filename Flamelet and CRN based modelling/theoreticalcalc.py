import cantera as ct
gas = ct.Solution('gri30.cti')
gas.TPY=300, ct.one_atm, {'CH4':1}
r1=ct.Reservoir(gas)
gas.TPY=673.15, ct.one_atm, {'O2':0.232,'N2':0.7547, 'CO2':0.00046, 'AR':0.0128}
r2= ct.Reservoir(gas)
gas.TPY=300, ct.one_atm, {'O2':0.232,'N2':0.7547, 'CO2':0.00046, 'AR':0.0128}
igr=ct.IdealGasReactor(gas)
igr.volume=
