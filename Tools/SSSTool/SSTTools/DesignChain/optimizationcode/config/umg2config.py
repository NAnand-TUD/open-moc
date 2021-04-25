from config import Config

class UMG2Config(Config):

    def __init__(self, name):
        Config.__init__(self)
        self.name = name
        self.set_default_properties()
        self.excludeconfig = ['log_filename', 'config_filename', 'name', 'output_filename']

    def set_default_properties(self):
        # Based on Mach 2.5 for original radius
        self.log_filename = 'umg2'+self.name+'.log'
        self.config_filename = 'umg2'+self.name+'.cfg'
        self.SPECS_FILE = 'stator_specs.out'
        self.FILE_NAME = 'stator'
        self.TE_thk = 0.0005
        self.RAD_FAC = 6
        self.TE_coords = '[1.0669938448821767,-0.5263474471065392]'
        self.Throat = 0.0012134982407772388
        self.MIN_NELEM = 1e5
        self.scale = 10
        self.BL_IN = 'FAC'
        self.BL_H = 0.005
        self.BL_FAC = 0.15