import datajoint as dj
from pipeline import map, experiment, shared

schema = dj.schema('pipeline_anatomy')

class Area (dj.Lookup):
     definition = """
          brain_area : varchar(12)   # short brain area name
     """
     contents = [
            ['V1'],
            ['P'],
            ['POR'],
            ['PM'],
            ['AM'],
            ['A'],
            ['RL'],
            ['AL'],
            ['LI'],
            ['LM']
     ]

class AreaMask (dj.Manual):
     definition = """
          # Area mask for each scan
          -> experiment.Scan
          -> Area
          -> shared.Field
          ---
          -> map.RetMap
          mask                     : mediumblob            # mask of area
     """

     # Convert MATLAB function here