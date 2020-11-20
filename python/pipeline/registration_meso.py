import datajoint as dj
from . import stack, anatomy, meso

schema = dj.schema('pipeline_meso')

@schema
class StackCoordinates(dj.Computed):
    definition = """ # centroids of each unit in motor/stack coordinate system

    -> meso.ScanSet          # animal_id, session, scan_idx, channel, field, segmentation_method, pipe_version
    -> stack.Registration.proj(session='scan_session')  # animal_id, stack_session, stack_idx, volume_id, session, scan_idx, field, stack_channel, scan_channel, registration_method
    """

    class UnitInfo(dj.Part):
        definition = """ # ScanSet.UnitInfo centroids mapped to stack coordinates

        -> master                       # this will add field and channels back
        -> meso.ScanSet.Unit
        ---
        stack_x         : float
        stack_y         : float
        stack_z         : float
        """

    def make(self, key):
        from scipy import ndimage

        # Get registration grid (px -> stack_coordinate)
        stack_key = {**key, 'scan_session': key['session']}
        field_res = (ScanInfo.Field & key).microns_per_pixel
        grid = (stack.Registration & stack_key).get_grid(type='affine',
                                                         desired_res=field_res)

        self.insert1(key)
        field_units = meso.ScanSet.UnitInfo & (meso.ScanSet.Unit & key)
        for unit_key, px_x, px_y in zip(*field_units.fetch('KEY', 'px_x', 'px_y')):
            px_coords = np.array([[px_y], [px_x]])
            unit_x, unit_y, unit_z = [ndimage.map_coordinates(grid[..., i], px_coords,
                                                              order=1)[0] for i in
                                      range(3)]
            StackCoordinates.UnitInfo.insert1({**key, **unit_key, 'stack_x': unit_x,
                                               'stack_y': unit_y, 'stack_z': unit_z})

@schema
class AreaMembership(dj.Computed):
    definition = """ # cell membership in visual areas according to stack registered retinotopy
    ret_hash             : varchar(32)                  # single attribute representation of retinotopic map key
    -> StackCoordinates
    ---
    """

    class UnitInfo(dj.Part):
        definition = """ # confidence in area assignment per unit according to stack coordinates
        -> master
        -> anatomy.Area
        -> StackCoordinates.UnitInfo
        ---
        confidence             : float                        # confidence in area assignment
        """
    @property
    def key_source(self):
        ret_rel = stack.Area.proj(ret_session='scan_session',
                                  ret_scan_idx='scan_idx',
                                  ret_channel='scan_channel')
        key_source = ret_rel * StackCoordinates
        heading_str = list(self.heading)

        return dj.U(*heading_str) & key_source

    def make(self,key):

        from scipy.interpolate import griddata

        mask_rel = stack.Area.Mask.proj('mask',
                                        ret_session='scan_session',
                                        ret_scan_idx='scan_idx',
                                        ret_channel='scan_channel')
        mask_keys, masks = (mask_rel & key).fetch('KEY', 'mask')

        fetch_str = ['x', 'y', 'um_width', 'um_height', 'px_width', 'px_height']
        stack_rel = stack.CorrectedStack.proj(*fetch_str, stack_session='session') & key
        cent_x, cent_y, um_w, um_h, px_w, px_h = stack_rel.fetch1(*fetch_str)

        stack_edges = np.array((cent_x - um_w / 2, cent_y - um_h / 2))
        stack_px_dims = np.array((px_w, px_h))
        stack_um_dims = np.array((um_w, um_h))

        stack_px_grid = np.meshgrid(*[np.arange(d) + 0.5 for d in stack_px_dims])

        ks, sxs, sys = (StackCoordinates.UnitInfo & key).fetch('KEY', 'stack_x', 'stack_y')
        pxs, pys = [np.array((coord - edge) * px_per_um) for coord, edge, px_per_um
                    in zip((sxs, sys), stack_edges, stack_px_dims / stack_um_dims)]

        unit_tups = []
        for mask_key, mask in zip(mask_keys, masks):
            grid_locs = np.array([grid.ravel() for grid in stack_px_grid]).T
            grid_vals = mask.ravel()
            grid_query = np.vstack((pxs, pys)).T

            confs = griddata(grid_locs, grid_vals, grid_query, method='nearest')
            mems = confs > 0

            unit_tups.append([{**mask_key, **k, 'confidence': conf} for k, conf in zip(np.array(ks)[mems], confs[mems])])

        self.insert1(key)
        self.UnitInfo.insert(np.concatenate(unit_tups),ignore_extra_fields=True)


@schema
class Func2StructMatching(dj.Computed):
    definition = """ # match functional masks to structural masks

    -> meso.ScanSet                  # animal_id, session, scan_idx, pipe_version, field, channel
    -> stack.FieldSegmentation.proj(session='scan_session') # animal_id, stack_session, stack_idx, volume_id, session, scan_idx, field, stack_channel, scan_channel, registration_method, stacksegm_channel, stacksegm_method
    ---
    key_hash        : varchar(32)       # single attribute representation of the key (used to avoid going over 16 attributes in the key)
    """

    class AllMatches(dj.Part):
        definition = """ # store all possible matches (one functional cell could match with more than one structural mask and viceversa)

        key_hash        : varchar(32)       # master key
        unit_id         : int               # functional unit id
        sunit_id        : int               # structural unit id
        ---
        iou             : float             # intersection-over-union of the 2-d masks
        """
        # Used key_hash because key using ScanSet.Unit, FieldSegmentation.StackUnit has
        # more than 16 attributes and MySQL complains. I added the foreign key constraints
        # manually

    class Match(dj.Part):
        definition = """ # match of a functional mask to a structural mask (1:1 relation)

        -> master
        -> meso.ScanSet.Unit
        ---
        -> stack.FieldSegmentation.StackUnit.proj(session='scan_session')
        iou             : float         # Intersection-over-Union of the 2-d masks
        distance2d      : float         # distance between centroid of 2-d masks
        distance3d      : float         # distance between functional centroid and structural centroid
        """

    def make(self, key):
        from .utils import registration
        from scipy import ndimage

        # Get caiman masks and resize them
        field_dims = (ScanInfo.Field & key).fetch1('um_height', 'um_width')
        masks = np.moveaxis((Segmentation & key).get_all_masks(), -1, 0)
        masks = np.stack([registration.resize(m, field_dims, desired_res=1) for m in
                          masks])
        scansetunit_keys = (meso.ScanSet.Unit & key).fetch('KEY', order_by='mask_id')

        # Binarize masks
        binary_masks = np.zeros(masks.shape, dtype=bool)
        for i, mask in enumerate(masks):
            ## Compute cumulative mass (similar to caiman)
            indices = np.unravel_index(np.flip(np.argsort(mask, axis=None), axis=0),
                                       mask.shape)  # max to min value in mask
            cumsum_mask = np.cumsum(mask[indices] ** 2) / np.sum(mask ** 2)  # + 1e-9)
            binary_masks[i][indices] = cumsum_mask < 0.9

        # Get structural segmentation and registration grid
        stack_key = {**key, 'scan_session': key['session']}
        segmented_field = (stack.FieldSegmentation & stack_key).fetch1('segm_field')
        grid = (stack.Registration & stack_key).get_grid(type='affine', desired_res=1)
        sunit_ids = (stack.FieldSegmentation.StackUnit & stack_key).fetch('sunit_id',
                                                                          order_by='sunit_id')

        # Create matrix with IOU values (rows for structural units, columns for functional units)
        ious = []
        for sunit_id in sunit_ids:
            binary_sunit = segmented_field == sunit_id
            intersection = np.logical_and(binary_masks, binary_sunit).sum(
                axis=(1, 2))  # num_masks
            union = np.logical_or(binary_masks, binary_sunit).sum(
                axis=(1, 2))  # num_masks
            ious.append(intersection / union)
        iou_matrix = np.stack(ious)

        # Save all possible matches / iou_matrix > 0
        self.insert1({**key, 'key_hash': hash_key_values(key)})
        for mask_idx, func_idx in zip(*np.nonzero(iou_matrix)):
            self.AllMatches.insert1({'key_hash': hash_key_values(key),
                                     'unit_id': scansetunit_keys[func_idx]['unit_id'],
                                     'sunit_id': sunit_ids[mask_idx],
                                     'iou': iou_matrix[mask_idx, func_idx]})

        # Iterate over matches (from best to worst), insert
        while iou_matrix.max() > 0:
            # Get next best
            best_mask, best_func = np.unravel_index(np.argmax(iou_matrix),
                                                    iou_matrix.shape)
            best_iou = iou_matrix[best_mask, best_func]

            # Get stack unit coordinates
            coords = (stack.FieldSegmentation.StackUnit & stack_key &
                      {'sunit_id': sunit_ids[best_mask]}).fetch1('sunit_z', 'sunit_y',
                                                                 'sunit_x', 'mask_z',
                                                                 'mask_y', 'mask_x')
            sunit_z, sunit_y, sunit_x, mask_z, mask_y, mask_x = coords

            # Compute distance to 2-d and 3-d mask
            px_y, px_x = ndimage.measurements.center_of_mass(binary_masks[best_func])
            px_coords = np.array([[px_y], [px_x]])
            func_x, func_y, func_z = [ndimage.map_coordinates(grid[..., i], px_coords,
                                                              order=1)[0] for i in
                                      range(3)]
            distance2d = np.sqrt((func_z - mask_z) ** 2 + (func_y - mask_y) ** 2 +
                                 (func_x - mask_x) ** 2)
            distance3d = np.sqrt((func_z - sunit_z) ** 2 + (func_y - sunit_y) ** 2 +
                                 (func_x - sunit_x) ** 2)

            self.Match.insert1({**key, **scansetunit_keys[best_func],
                                'sunit_id': sunit_ids[best_mask], 'iou': best_iou,
                                'distance2d': distance2d, 'distance3d': distance3d})

            # Deactivate match
            iou_matrix[best_mask, :] = 0
            iou_matrix[:, best_func] = 0

