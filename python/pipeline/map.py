import datajoint as dj
from . import mouse

schema = dj.schema('pipeline_map')

@schema
class RetMap (dj.Manual):
     definition = """
          # Retinotopy map
          -> mouse.Mouse
          ret_idx : smallint        # retinotopy map index for each animal
          ---
     """

     # Convert MATLAB function here