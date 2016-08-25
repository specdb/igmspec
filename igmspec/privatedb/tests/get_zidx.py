from linetools import utils as ltu
coords = [ltu.radec_to_coord(jc) for jc in ['J001605.89+005654.3','J001607.27+005653.1']]
from igmspec.igmspec import IgmSpec
from igmspec import cat_utils as icu
igmsp = IgmSpec()
# PUT A PDB IN cat_utils
zem = icu.zem_from_radec([icoord.ra.value for icoord in coords],
                         [icoord.dec.value for icoord in coords], igmsp.idb.hdf)
# 161121, 161130
