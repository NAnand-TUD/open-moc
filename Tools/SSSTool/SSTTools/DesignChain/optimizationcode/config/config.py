class Config(object):
    """
    The config objects are the python representation of the configuration files for each executable
    """
    def __init__(self):
        """
        :ivar str output_filename: initial value: ''
        :ivar str config_filename: initial value: ''
        :ivar str log_filename: initial value: ''
        :ivar  excludeconfig: initial value: []
        :vartype excludeconfig: list[str]
        """
        self.output_filename = ""
        self.config_filename = ""
        self.log_filename = ""
        self.excludeconfig = []
        
    def get_config_filecontent(self):
        """
        Convert the object attributes (except the attributes given in excludeconfig) into text ready to be written into configuration file
        
        :return: filecontent
        :rtype: str
        """
        lines =[]
        self.excludeconfig.append("excludeconfig")
        orderedlist = sorted(self.__dict__)
        for key in orderedlist	:
            if not (key in self.excludeconfig):
                lines.append(key + "=" + str(self.__dict__[key]))
        return "\n".join(lines)
