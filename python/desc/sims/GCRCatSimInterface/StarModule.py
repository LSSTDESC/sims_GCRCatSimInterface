import numpy as np
from sqlalchemy import text
from lsst.sims.utils import htmModule as htm
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catalogs.db import ChunkIterator
from lsst.sims.utils import _angularSeparation

__all__ = ["DC2StarObj"]


class DC2StarObj(CatalogDBObject):

    objid = 'dc2_stars'
    tableid = 'stars'
    objectTypeId = 4
    idColKey = 'simobjid'
    raColKey = 'ra'
    decColKey = 'decl'

    columns = [('id', 'simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('galacticAv', '3.1*ebv'),
               ('variabilityParameters', 'varParamStr', str, 256)]

    def _final_pass(self, chunk):
        """
        Only keep objects that are inside the field of view
        """
        if len(chunk) == 0:
            return chunk

        dd = _angularSeparation(self._obs_metadata._pointingRA,
                                self._obs_metadata._pointingDec,
                                chunk['raJ2000'], chunk['decJ2000'])

        valid = np.where(dd<=self._obs_metadata._boundLength)
        return chunk[valid]

    def query_columns(self, colnames=None, chunk_size=None,
                      obs_metadata=None, constraint=None, limit=None):

        query = self._get_column_query(colnames)
        self._obs_metadata = obs_metadata

        half_space = htm.halfSpaceFromRaDec(obs_metadata.pointingRA,
                                            obs_metadata.pointingDec,
                                            obs_metadata.boundLength)

        htmid_bounds = half_space.findAllTrixels(6)
        htmid_clause = '('
        for bound in htmid_bounds:
            if htmid_clause != '(':
                htmid_clause += ' OR'
            if bound[0] == bound[1]:
                htmid_clause += '(htmid_6 == %d)' % bound[0]
            else:
                htmid_clause += '(htmid_6 <= %d AND htmid_6 >= %d)'  % (bound[1], bound[0])

        htmid_clause += ')'

        query = query.filter(text(htmid_clause))

        return ChunkIterator(self, query, chunk_size)

