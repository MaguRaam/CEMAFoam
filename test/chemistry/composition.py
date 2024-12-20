import cantera as ct

gas = ct.Solution("gri30.yaml")
air = "O2:0.21,N2:0.79"
fuel = "CH4:1"
gas.TPX = 1000, 1e6, air
gas()

gas.set_equivalence_ratio(phi=1, fuel=fuel, oxidizer=air, basis="mole")
gas()