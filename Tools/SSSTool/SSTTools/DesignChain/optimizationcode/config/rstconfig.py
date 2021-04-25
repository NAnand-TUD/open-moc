from config import Config

class RSTConfig(Config):

    def __init__(self, name):
        Config.__init__(self)
        self.name = name
        self.set_default_properties()
        self.excludeconfig = ['log_filename', 'config_filename', 'name', 'output_filename', 'AreaRatio']

    def set_default_properties(self):
        self.log_filename = 'rst'+self.name+'.log'
        self.config_filename = 'rst'+self.name+'.cfg'
        self.nBlades = 18
        self.flowAngle = 80.4 
        self.radiusOut = 0.11908
        self.outFileMoc = 'nozzle_coords.out'
        self.kernelRadiusRatio = 10.0
        self.ActualRin = 0.1685
        self.ActualRout = 0.11908
        self.RealRout = 0.11
        self.ScaleNoz = 1
        self.TEminT = 0.0005
        self.ScaleMax = 50
        #self.Throat = 1e-3
        self.AreaRatio = 25
        self.Mode = 'Sim' 
        self.OptiDiv = 'Yes'
        self.ScaleMesh = 10
        self.nPointsScale = 2
        self.CoordsName =  'stator_coords.out'
        self.UMG2Name = 'stator'
        self.SpecsName = 'stator_specs.out'



