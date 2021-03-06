"""
Storage.py

Created by Dave Williams on 2009-06-28.
"""

import yaml
import cPickle as pickle
import datetime
import warnings

class Storage():
    """Interface with a stored set of crossbridge data.
        
    Use list() to see what properties are currently stored
        get(xb_property) to access properties
        write(xb_property, new_value) change properties
        save() to save to disk"""
    
    def __init__(self, xbtype, config=None, x_range=None, y_range=None, protocol="pickle"):
        """ Parse the variables and deal with the four cases in which 
        Storage will be used, here is an outline of the process:
         - construct the file name from the number of springs and protocol
         - try to read in the file
         - if that fails (because the file doesn't exist)
            - and params were passed (so we will be wanting to write data 
              to disk later), then create a new file with the passed params
            - and no params were passed (because we are really only reading 
              the file), then throw an exception and get out of dodge
         - if params were passed and any of them don't match our file, trash
           the old ones from the file and use the new ones
        """
        self.protocol = protocol
        self.file_name = str(xbtype)+"spring."
        ## Read in file, or create if needed
        try:
            if self.protocol == "pickle":
                stream = open(self.file_name+"pkl", 'r')
                self.data = pickle.load(stream)
            elif self.protocol == "yaml":
                stream = open(self.file_name+"yml", 'r')
                self.data = yaml.load(stream)
            else:
                warnings.warn("Unrecognized file protocol, aborting")
                return
            stream.close()
        except IOError: # File doesn't exist
            if config is not None: # Passed params, in write mode
                msg = ("\n Storage file " + self.file_name + 
                " not found, creating one with passed parameters")
                warnings.warn(msg)
                self.data = {
                    "xbtype": xbtype,
                    "config": config,
                    "x_range": x_range,
                    "y_range": y_range
                }
                self.save()
            else: # No passed params, in read mode
                msg = ("\n Storage file " + self.file_name + 
                " not found, can't read in data to continue")
                raise Exception(msg)
        # In case some key parameter has changed, trash the old file
        if (config is not None) and (self.get('config') != config):
            msg = "\n Config has changed, trashing old data and starting anew"
            warnings.warn(msg)
            self.__trash__(xbtype, config, x_range, y_range)
        if (config is not None) and (self.get('x_range') != x_range or 
            self.get('y_range') != y_range):
            msg = "\n Range has changed, trashing old data and starting anew"
            warnings.warn(msg)
            self.__trash__(xbtype, config, x_range, y_range)
    
    def __trash__(self, xbtype, config, x_range, y_range):
        """Trash the stored data"""
        for key in self.data.keys():
            del self.data[key]
        self.data = {
            "xbtype": xbtype,
            "config": config,
            "x_range": x_range,
            "y_range": y_range
            }
        
    
    def list(self):
        """List the current xb_properties"""
        return self.data.keys()
    
    def get(self, xb_property):
        """Return the stored property for the file, if it exists."""
        if self.data.has_key(xb_property):
            return self.data.get(xb_property)
        else:
            raise Exception("Storage was asked for an unknown xb property")
    
    def write(self, xb_property, new_value):
        """Write a property to the current data; doesn't save to disk."""
        if xb_property not in ("xbtype", "config", 
                                "x_range", "y_range", "timestamp"):
            self.data[xb_property] = new_value
        else:
            warnings.warn("Don't set base attributes, should be instantiated")
    
    def save(self, protocol=None):
        """Save the current data to a file on disk."""
        if protocol is None:
            protocol = self.protocol # Dance with the one who brought ya
        self.data['timestamp'] = datetime.datetime.today()
        if protocol == "pickle":
            stream = open(self.file_name+"pkl", 'w')
            pickle.dump(self.data, stream)
        elif protocol == "yaml":
            stream = open(self.file_name+"yml", 'w')
            stream.write(yaml.dump(self.data, indent=4))
        else:
            warnings.warn("Unrecognized save format, data not saved")
        stream.close()
    

