from lsst.sims.utils import htmModule

__all__ = ["DC2_bounds"]

class DC2_bounds(object):

    def __init__(self):
        ne_corner = (71.46, -27.25)
        nw_corner = (52.25, -27.25)
        se_corner = (73.79, -44.33)
        sw_corner = (49.92, -44.33)
        pt_inside = (60.0, -30.0)

        self._hs_list = []
        self._hs_list.append(htmModule.halfSpaceFromPoints(ne_corner,
                                                           nw_corner,
                                                           pt_inside))

        self._hs_list.append(htmModule.halfSpaceFromPoints(ne_corner,
                                                           se_corner,
                                                           pt_inside))

        self._hs_list.append(htmModule.halfSpaceFromPoints(se_corner,
                                                           sw_corner,
                                                           pt_inside))

        self._hs_list.append(htmModule.halfSpaceFromPoints(sw_corner,
                                                           nw_corner,
                                                           pt_inside))

    @property
    def hs_list(self):
        return self._hs_list
