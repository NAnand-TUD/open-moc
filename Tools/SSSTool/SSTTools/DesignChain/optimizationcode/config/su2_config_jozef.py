from config import Config

class SU2Config(Config):

    def __init__(self, name):
        Config.__init__(self)
        self.name = name
        self.set_default_properties()
        self.excludeconfig = ['log_filename', 'config_filename', 'name', 'output_filename']

    def set_default_properties(self):
        self.log_filename = 'su2_per' + self.name + '.log'
        self.config_filename = 'su2_per' + self.name + '.cfg'
        self.PHYSICAL_PROBLEM = 'RANS'
        self.KIND_TURB_MODEL = 'SST'
        self.MATH_PROBLEM = 'DIRECT'
        self.RESTART_SOL = 'NO'
        self.MACH_NUMBER = 0.05
        self.AoA = 0.0
        self.FREESTREAM_PRESSURE = 3195000.0
        self.FREESTREAM_TEMPERATURE = 587.65
        self.FREESTREAM_DENSITY = 0.000488#98.9258 #1.2886
        self.FREESTREAM_TURBULENCEINTENSITY = 0.1
        self.FREESTREAM_TURB2LAMVISCRATIO = 100.0
        self.FREESTREAM_OPTION = 'TEMPERATURE_FS'
        self.INIT_OPTION = 'TD_CONDITIONS'
        self.REF_DIMENSIONALIZATION = 'DIMENSIONAL' #'FREESTREAM_PRESS_EQ_ONE'
        self.FLUID_MODEL = 'PR_GAS'
        self.GAMMA_VALUE = 1.06 #1.35 #1.06
        self.GAS_CONSTANT = 90.23
        self.CRITICAL_TEMPERATURE = 591.75
        self.CRITICAL_PRESSURE = 4126300.0
        self.ACENTRIC_FACTOR = 0.2657
        self.VISCOSITY_MODEL = 'CONSTANT_VISCOSITY'
        self.MU_CONSTANT = 6.799E-06#1.519E-05 #1.644E-4
        self.MU_REF = 1.716E-5
        self.MU_T_REF = 273.15
        self.SUTHERLAND_CONSTANT = 110.4
        self.CONDUCTIVITY_MODEL = 'CONSTANT_CONDUCTIVITY'
        self.KT_CONSTANT = 0.01014#0.04602 #0.08029
        self.MARKER_HEATFLUX = '(wall1, 0.0)'
        self.MARKER_GILES = '( inflow, TOTAL_CONDITIONS_PT, 3195000.0, 587.65, -1.0, 0.0, 0.0, 0.0, 0.0, outflow, STATIC_PRESSURE, 0.8E+05, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)'
        self.MARKER_RIEMANN = '( NONE )'
        self.AVERAGE_PROCESS_KIND = 'MIXEDOUT'
        self.TURBOMACHINERY_KIND = 'CENTRIPETAL'
        self.NUM_SPANWISE_SECTIONS = 1
        self.RAMP_OUTLET_PRESSURE = 'YES'
        self.RAMP_OUTLET_PRESSURE_COEFF = (2000000.0, 10.0, 1000)
        self.MARKER_PERIODIC = '(periodic1, periodic2, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0)'
        self.MARKER_PLOTTING = '(outflow)'
        self.MARKER_MONITORING = '(wall1)'
        self.MARKER_TURBOMACHINERY = '(inflow, outflow)'
        self.KIND_ADAPT = 'PERIODIC'
        self.GRID_MOVEMENT = 'NO'
        self.GRID_MOVEMENT_KIND = 'ROTATING_FRAME'
        self.MACH_MOTION = 0.35
        self.ROTATION_RATE_X = '0.0 0.0'
        self.ROTATION_RATE_Y = '0.0 0.0'
        self.ROTATION_RATE_Z = '0.0 1990.0'
        self.NUM_METHOD_GRAD = 'WEIGHTED_LEAST_SQUARES'
        self.CFL_NUMBER = 40.0
        self.CFL_ADAPT = 'NO'
        self.CFL_ADAPT_PARAM = (0.3, 0.5, 1.0, 1000.0)
        self.LINEAR_SOLVER = 'FGMRES'
        self.LINEAR_SOLVER_PREC = 'LU_SGS'
        self.LINEAR_SOLVER_ERROR = 1E-4
        self.LINEAR_SOLVER_ITER = 5
        self.MGLEVEL = 0
        self.MG_PRE_SMOOTH = (1, 2, 3, 3)
        self.MG_POST_SMOOTH = (0, 0, 0, 0)
        self.MG_CORRECTION_SMOOTH = (0, 0, 0, 0)
        self.MG_DAMP_RESTRICTION = 0.75
        self.MG_DAMP_PROLONGATION = 0.75
        self.REF_ELEM_LENGTH = 0.1
        self.LIMITER_COEFF = 0.5
        self.LIMITER_ITER = 999999
        self.CONV_NUM_METHOD_FLOW = 'ROE'
        self.SPATIAL_ORDER_FLOW = '2ND_ORDER_LIMITER'
        self.ENTROPY_FIX_COEFF = 0.001
        self.AD_COEFF_FLOW = (0.15, 0.5, 0.02)
        self.SLOPE_LIMITER_FLOW = 'VAN_ALBADA'
        self.TIME_DISCRE_FLOW = 'EULER_IMPLICIT'
        self.RELAXATION_FACTOR_FLOW = 1.0
        self.CONV_NUM_METHOD_TURB = 'SCALAR_UPWIND'
        self.SPATIAL_ORDER_TURB = '1ST_ORDER'
        self.SLOPE_LIMITER_TURB = 'VENKATAKRISHNAN'
        self.TIME_DISCRE_TURB = 'EULER_IMPLICIT'
        self.CFL_REDUCTION_TURB = 0.01
        self.RELAXATION_FACTOR_TURB = 1.0
        self.DV_KIND = 'SCALE'
        self.DV_MARKER = '(airfoil, wall1, periodic1, periodic2, inflow, outflow)'
        self.DV_PARAM = '1, 1'
        self.DV_VALUE = 0.001
        self.EXT_ITER = 10001
        self.CONV_CRITERIA = 'RESIDUAL'
        self.RESIDUAL_FUNC_FLOW = 'RHO_ENERGY'
        self.RESIDUAL_REDUCTION = 6
        self.RESIDUAL_MINVAL = -16
        self.STARTCONV_ITER = 10
        self.CAUCHY_ELEMS = 100
        self.CAUCHY_EPS = 1E-6
        self.CAUCHY_FUNC_FLOW = 'DRAG'
        self.MESH_FILENAME = 'su2mesh.su2'
        self.MESH_FORMAT = 'SU2'
        self.MESH_OUT_FILENAME = 'su2mesh_per.su2'
        self.SOLUTION_FLOW_FILENAME = 'restart_flow.dat'
        self.CONV_FILENAME = 'history'
        self.RESTART_FLOW_FILENAME = 'restart_flow.dat'
        self.VOLUME_FLOW_FILENAME = 'flow'
        self.SURFACE_FLOW_FILENAME = 'surface_flow'
        self.WRT_SOL_FREQ = 50
        self.WRT_CON_FREQ = 1