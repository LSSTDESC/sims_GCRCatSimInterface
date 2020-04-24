import os
import glob
from lsst.afw import cameraGeom
from lsst.obs.lsstSim import LsstSimMapper
from desc.imsim import InstCatTrimmer

def trim_file_name(sensor_id, outdir='.'):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    detname = 'R{}{}_S{}{}'.format(*[_ for _ in sensor_id if _.isdigit()])
    return os.path.join(outdir, f'{detname}_instcat.txt')

def write_instcat(commands, obj_lines, outfile):
    with open(outfile, 'w') as output:
        for command in commands:
            output.write("{}\n".format(command.strip()))
        for line in obj_lines:
            output.write("{}\n".format(line.strip()))

camera = LsstSimMapper().camera
sensor_list = [det.getName() for det in camera
               if det.getType() == cameraGeom.SCIENCE]

instcats = sorted(glob.glob('00*/phosim_cat*.txt'))
for instcat in instcats:
    folder = os.path.dirname(instcat)
    outdir = f'instcats_{folder}'
    os.makedirs(outdir, exist_ok=True)
    trimmer = InstCatTrimmer(instcat, sensor_list, log_level='DEBUG')
    for sensor_id in sensor_list:
        outfile = trim_file_name(sensor_id, outdir)
        if os.path.isfile(outfile):
            continue
        write_instcat(trimmer.command_lines, trimmer[sensor_id], outfile)
