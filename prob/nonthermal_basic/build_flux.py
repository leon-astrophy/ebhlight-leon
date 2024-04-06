################################################################################
#                                                                              #
#  NONTHERMAL TESTING                                                          #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'nonthermal_basic'


### COMPILE TIME PARAMETERS ###
# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 256)
bhl.config.set_cparm('N2TOT', 1)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', False)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# ELECTRONS
bhl.config.set_cparm('ELECTRONS', True)
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', False)
bhl.config.set_cparm('BETA_HEAT', True)
bhl.config.set_cparm('COULOMB', True)

# NONTHERMAL
bhl.config.set_cparm('NONTHERMAL', True)
bhl.config.set_cparm('SYNCHROTRON', True)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('ADIABTIC_SCALING', True)
bhl.config.set_cparm('SKIP_ADIAB', True)
bhl.config.set_cparm('SKIP_VISCOUS', True)
bhl.config.set_cparm('SKIP_COOLING', True)
bhl.config.set_cparm('LOGDUMPING', True)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)


### RUNTIME PARAMETERS ###
bhl.config.set_rparm('tf', 'double', default = 0.25)
bhl.config.set_rparm('dt', 'double', default = 0.25e-6)
bhl.config.set_rparm('gam', 'double', default = 1.4)
bhl.config.set_rparm('L_unit', 'double', default = 1.5032e6)
bhl.config.set_rparm('M_unit', 'double', default = 1.42034e12)
bhl.config.set_rparm('DTd', 'double', default = 0.25e-1)
bhl.config.set_rparm('DTl', 'double', default = 0.25e-2)
bhl.config.set_rparm('DTr', 'double', default = 12837612)
bhl.config.set_rparm('test_type', 'integer', default = 1)

bhl.config.set_cparm('PLAW', 3.5)
bhl.config.set_rparm('plaw', 'double', default = 3.5)


### CONFIGURE AND COMPILE  ###
bhl.build(PROB)

