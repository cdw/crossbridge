"""
Storage.py

Created by Dave Williams on 2009-06-28.
"""

import yaml
import datetime
import warnings

class Storage():
    """Interface with a stored set of crossbridge data.
        
    Use list() to see what properties are currently stored
        get(xb_property) to access properties
        write(xb_property, new_value) change properties
        save() to save to disk"""
    
    def __init__(self, xbtype, config, x_range, y_range):
        self.xbtype = xbtype
        self.config = config
        self.x_range = x_range
        self.y_range = y_range
        self.base_types = set(['xbtype', 'config', 'x_range', 'y_range'])
        self.data_file_name = str(self.xbtype)+"spring.yaml"
        try:
            stream = open(self.data_file_name, 'r')
            self.data = yaml.load(stream)
            stream.close()
        except IOError:
            warnings.warn("File "+self.data_file_name+" not found, creating")
            stream = open(self.data_file_name,"w")
            self.data = yaml.load("""config: None""")
            stream.close()
            return
        ## Check stored data, make sure we aren't changing too much
        if self.get('config') != self.config:
            warnings.warn("Using a different XB config than before.\n \
            Gonna trash the old data and only use what is added next.")
            self.__trash__()
        if (self.get('x_range') != self.x_range or 
            self.get('y_range') != self.y_range):
            warnings.warn("Using a different range than before.\n \
            Gonna trash the old data and only use what is added next.")
            self.__trash__()
    
    def __trash__(self):
        """Trash the stored data"""
        for key in self.data.keys():
            del self.data[key]
        self.data['config'] = self.config
        self.data['xbtype'] = self.xbtype
        self.data['x_range'] = self.x_range
        self.data['y_range'] = self.y_range
        self.data['modified'] = datetime.datetime.today()
    
    def list(self):
        """List the current xb_properties"""
        return self.data.keys()
    
    def get(self, xb_property):
        """Return the stored property for the file, if it exists."""
        if self.data.has_key(xb_property):
            return self.data.get(xb_property)
        else:
            return -1
    
    def write(self, xb_property, new_value):
        """Write a property to the current data; doesn't save to disk."""
        if xb_property not in self.base_types:
            self.data[xb_property] = new_value
        else:
            warnings.warn("Don't set base attributes, must be instantiated")
    
    def save(self):
        """Save the current data to a YAML file on disk."""
        stream = open(self.data_file_name, 'w')
        stream.write(yaml.dump(self.data))
        stream.close()
    
