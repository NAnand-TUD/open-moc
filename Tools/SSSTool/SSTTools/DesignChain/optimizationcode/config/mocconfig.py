from config import Config

class MOCConfig(Config):

    def __init__(self, name):
        Config.__init__(self)
        self.name = name
        self.set_default_properties()
        self.excludeconfig = ['log_filename', 'config_filename', 'name', 'output_filename', 'File_Name_Plot']

    def set_default_properties(self):
        self.log_filename = 'moc'+self.name+'.log'
        self.config_filename = 'moc'+self.name+'.cfg'
        self.NozzleType = "PLANAR"
        self.GAS_EQU = "CoolProp"
        self.FLDNAME = "Toluene"
        self.THROAT_PROP = "ITER"
        self.VrhoGuess = 0.03
        self.TGuess = 580
        self.PGuess = 20e5
        self.Tc = 591.75 #591.79
        self.Pc = 4126300.0 #2.09e6
        self.Vc = 0.000316
        self.M = 0.09214
        self.Ho = 676.019
        self.so = 0.9
        self.gamma = 1.055
        self.To = 580
        self.Po = 32e5#20.4e5
        self.R = 8.314
        self.rho_t = 10
        self.y_t = 1
        self.rho_d = 10
        self.n = 10
        self.tolerance_v = 1e-6
        self.tolerance_x = 1e-6
        self.dtau = 0.5
        self.n_ref = 40
        self.Noz_Design_Mach = 2.5
        self.File_Name_NProp = "nozzle_prop.out"
        self.File_Name_NCods = "nozzle_coords.out"
        self.File_Name_Plot = "characteristics.png"




