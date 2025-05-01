PSUEDO CODE:
============================================
1/ init chem system:
   - load mech (e.g., gri-mech 3.0)
   - set init state: temp, pressure, species comp

2/ define ode system:
   - num eq = num species + 1 (temp)
   - for each step:
     1/ update thermo state (curr temp, mass fracs)
     2/ compute net prod rates
     3/ calc derivs:
        - dY_i/dt = (omega_i * mw_i) / density
        - dT/dt = -hdot / (density * cp)

3/ set init conds:
   - mass fracs, temp

4/ init cvode:
   - set tolerances
   - give ode + init conds

5/ integrate:
   - solve odes t=0 to t_end

6/ get + output results:
   - extract solution (mass fracs, temp)
   - print final temp

7/ handle errors, cleanup:
   - catch, show errors
   - delete integrator, free resources



CODE LOGIC:
============================================
class ReactorODEs:
    ctor(sol):
        - init thermo + kinetic models from sol
        - store p, num species, alloc state vec

    neq():
        - return num eq (species + temp)

    eval(t, y, ydot, p):
        - update thermo state from y (mass fracs, temp)
        - compute net prod rates (omega)
        - calc derivs:
            - dY_i/dt = (omega_i * mw_i) / density
            - dT/dt = -hdot / (density * cp)

    getInitialConditions(t0, y):
        - set init mass fracs, temp in y

    getState(y):
        - get curr state (mass fracs, temp)
        
        
main:
    1/ init chem system:
        - load gri30.yaml
        - set T=1000K, P=1atm, comp="H2:2,O2:1,N2:4"

    2/ make ode system:
        - new ReactorODEs(chem sys)

    3/ setup integrator:
        - new cvode integrator
        - init t=0, ode system
        - set tolerances

    4/ solve odes:
        - integrate t=0 to t=1.0s

    5/ get + output results:
        - extract solution (mass fracs, temp)

