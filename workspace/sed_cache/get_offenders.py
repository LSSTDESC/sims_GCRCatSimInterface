import h5py
import os
import numpy as np

def get_offenders(av, rv):
    small_rv = rv<0.1
    med_rv_big_av = np.logical_and(rv<1.0,av>1.0)
    small_av = av<0.0
    return np.logical_or(small_rv,
                    np.logical_or(med_rv_big_av, small_av))

in_dir = os.path.join(os.environ['SCRATCH'], 'sed_cache')
assert os.path.isdir(in_dir)

disk_name = os.path.join(in_dir, 'disk_10451.h5')
assert os.path.isfile(disk_name)

bulge_name = os.path.join(in_dir, 'bulge_10451.h5')
assert os.path.isfile(bulge_name)

disk = h5py.File(disk_name, 'r')
bulge = h5py.File(bulge_name, 'r')

disk_gid = disk['galaxy_id'].value
bulge_gid = bulge['galaxy_id'].value

np.testing.assert_array_equal(disk_gid, bulge_gid)

disk_av = disk['A_v'].value
disk_rv = disk['R_v'].value
bulge_av = bulge['A_v'].value
bulge_rv = bulge['R_v'].value

offensive_disks = get_offenders(disk_av, disk_rv)
offensive_bulges = get_offenders(bulge_av, bulge_rv)

offenders = np.where(np.logical_or(offensive_disks, offensive_bulges))
print('%d offenders' % len(offenders[0]))

out_disk = os.path.join(in_dir, 'disk_offenders_10451.h5')
assert not os.path.isfile(out_disk)
d = h5py.File(out_disk, 'w')
for kk in disk.keys():
    d.create_dataset(kk, data=disk[kk].value[offenders])
d.close()

out_bulge = os.path.join(in_dir, 'bulge_offenders_10451.h5')
assert not os.path.isfile(out_bulge)
b = h5py.File(out_bulge, 'w')
for kk in bulge.keys():
    b.create_dataset(kk, data=bulge[kk].value[offenders])
b.close()
