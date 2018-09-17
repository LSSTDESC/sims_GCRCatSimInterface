# This script can go through an instance catalog and apply some post-process
# corrections in order to avoid having to regenerate the whole catalog
import os
import numpy as np
import argparse
import contextlib
import gzip

@contextlib.contextmanager
def fopen(filename, **kwds):
    """
    Return a file descriptor-like object that closes the underlying
    file descriptor when used with the with-statement.
    Parameters
    ----------
    filename: str
        Filename of the instance catalog.
    **kwds: dict
        Keyword arguments to pass to the gzip.open or open functions.
    Returns
    -------
    generator: file descriptor-like generator object that can be iterated
        over to return the lines in a file.
    """
    abspath = os.path.split(os.path.abspath(filename))[0]
    try:
        if filename.endswith('.gz'):
            fd = gzip.open(filename, **kwds)
        else:
            fd = open(filename, **kwds)
        yield fopen_generator(fd, abspath, **kwds)
    finally:
        fd.close()

def fopen_generator(fd, abspath, **kwds):
    """
    Return a generator for the provided file descriptor that knows how
    to recursively read in instance catalogs specified by the
    includeobj directive.
    """
    with fd as input_:
        for line in input_:
            if not line.startswith('includeobj'):
                yield line
            else:
                filename = os.path.join(abspath, line.strip().split()[-1])
                with fopen(filename, **kwds) as my_input:
                    for line in my_input:
                        yield line


def metadata_from_file(file_name):
    """
    Read in the InstanceCatalog specified by file_name.
    Return a dict of the header values from that InstanceCatalog.
    Simpler version than the ImSim code
    """
    input_params = {}
    with fopen(file_name, mode='rt') as in_file:
        for line in in_file:
            if line[0] == '#':
                continue

            params = line.strip().split()

            if params[0] == 'object':
                break

            float_val = float(params[1])
            int_val = int(float_val)
            if np.abs(float_val-int_val)>1.0e-10:
                val = float_val
            else:
                val = int_val
            input_params[params[0]] = val

    commands = dict(((key, value) for key, value in input_params.items()))
    return commands



def fix_disk_knots(in_instcat_disk, in_instcat_knots,
         out_instcat_disk, out_instcat_knots):

    # Use .fopen to read in the command and object lines from the
    # instance catalog.
    count_extinction = 0
    with fopen(in_instcat_disk, mode='rt') as input_disk,   \
         fopen(in_instcat_knots, mode='rt') as input_knots,  \
         open(out_instcat_disk, 'w') as output_disk, \
         open(out_instcat_knots, 'w') as output_knots:

        # We first go through the knots catalog because some entries are missing
        # compared to the full bulge/disk catalog (faint knots have already been
        # removed at the instance catalog creation level)
        for line_knots in input_knots:
            # Extract the galaxy ID for that knots component
            tokens_knots = line_knots.strip().split()
            id_knots = int(tokens_knots[1]) >> 10

            found = False
            # Loop through the disk catalogs
            for line_disk in input_disk:
                tokens_disk = line_disk.strip().split()
                id_disk = int(tokens_disk[1]) >> 10

                # Extracting the extinction parameters:
                if tokens_disk[17].lower() != 'none':
                    disk_internal_av = float(tokens_disk[18])
                    disk_internal_rv = float(tokens_disk[19])
                else:
                    disk_internal_av = 0
                    disk_internal_rv = 0

                # If the galaxy is offensive, clip the av and rv values
                if disk_internal_av < 0 or disk_internal_rv < 1:
                    #print('Fixing offensive disk %d with av/rv extinction: %f/%f'%(id_disk,disk_internal_av, disk_internal_rv))
                    disk_internal_av = np.clip(disk_internal_av,0.0,None)
                    disk_internal_rv = np.clip(disk_internal_rv,1.0,None)
                    tokens_disk[18] = ("%.7f"%disk_internal_av).rstrip('0')
                    tokens_disk[19] = ("%.7f"%disk_internal_rv).rstrip('0')
                    count_extinction += 1

                if id_disk == id_knots:
                    found=True
                    break
                else:
                    line_disk = ' '.join(tokens_disk)
                    output_disk.write(line_disk.strip()+'\n')

            if not found:
                print("ERROR: object ids do not match between input knots and disks catalogs")
                exit(-1)

            # Get total flux
            magnorm_disk = np.float(tokens_disk[4])
            magnorm_knots = np.float(tokens_knots[4])
            total_flux = 10.**(-magnorm_disk/2.5) + 10.**(-magnorm_knots/2.5)
            knots_flux_ratio = 10.**(-magnorm_knots/2.5) / total_flux

            # Apply smooth flux ration scaling
            size = np.float(tokens_disk[13])
            knots_flux_ratio = knots_flux_ratio * 0.5*(1 - np.tanh(np.log(size/1.5)))

            magnorm_disk = -2.5*np.log10((1-knots_flux_ratio)*total_flux)
            magnorm_knots = -2.5*np.log10(knots_flux_ratio*total_flux)

            # Update the entry
            tokens_disk[4] = ("%.9f"%magnorm_disk).rstrip('0')
            tokens_knots[4] = ("%.9f"%magnorm_knots).rstrip('0')
            # Making sure that the extinction paramters remain the same between disk and knots
            tokens_knots[18] = tokens_disk[18]
            tokens_knots[19] = tokens_disk[19]
            line_disk = ' '.join(tokens_disk)
            line_knots = ' '.join(tokens_knots)

            # Write the catalogs
            output_disk.write(line_disk.strip()+'\n')
            output_knots.write(line_knots.strip()+'\n')

    print("Fixed extinction for %d disks"%count_extinction)


def fix_bulge(in_instcat_bulge, out_instcat_bulge):
    count_extinction = 0
    with fopen(in_instcat_bulge, mode='rt') as input_bulge,  \
         open(out_instcat_bulge, 'w') as output_bulge:

        # Now we go through the bulge catalog independently, fixing the extinction
        for line_bulge in input_bulge:
            tokens_bulge = line_bulge.strip().split()
            id_bulge = int(tokens_bulge[1]) >> 10

            # Extracting the extinction parameters:
            if tokens_bulge[17].lower() != 'none':
                internal_av = float(tokens_bulge[18])
                internal_rv = float(tokens_bulge[19])
            else:
                internal_av = 0
                internal_rv = 0

            # If the galaxy is offensive, clip the av and rv values
            if internal_av < 0 or internal_rv < 1:
                #print('Fixing offensive bulge %d with av/rv extinction: %f/%f'%(id_bulge,internal_av, internal_rv))
                internal_av = np.clip(internal_av,0.0,None)
                internal_rv = np.clip(internal_rv,1.0,None)
                tokens_bulge[18] = ("%.9f"%internal_av).rstrip('0')
                tokens_bulge[19] = ("%.9f"%internal_rv).rstrip('0')
                count_extinction += 1

            line_bulge = ' '.join(tokens_bulge)
            output_bulge.write(line_bulge.strip()+'\n')
    print("Fixed extinction for %d bulge"%count_extinction)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Instance catalog crawler applying corrections in post-processing')
    parser.add_argument('input_cat', type=str,help='Phosim instance catalog to process')
    parser.add_argument('output_path', type=str, help='Directory in which to store the corrected catalog')
    args = parser.parse_args()

    # Find the visit id
    metadata = metadata_from_file(args.input_cat)
    visitID = metadata['obshistid']
    input_path = args.input_cat.split('/')[:-1]
    input_path = '/'.join(input_path)
    print('Copying catalog from %s'%input_path)

    # Create output directory
    output_path= os.path.join(args.output_path,'{0:08d}'.format(int(visitID)))
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Copy over the content of the instance catalog
    os.system("cp -arv %s/* %s"%(input_path, output_path))

    # Processes catalogs
    input_disk=output_path+'/disk_gal_cat_%d.txt.gz'%visitID
    input_bulge=output_path+'/bulge_gal_cat_%d.txt.gz'%visitID
    input_knots=output_path+'/knots_cat_%d.txt.gz'%visitID

    output_disk=output_path+'/disk_gal_cat_%d.txt'%visitID
    output_bulge=output_path+'/bulge_gal_cat_%d.txt'%visitID
    output_knots=output_path+'/knots_cat_%d.txt'%visitID

    print('Processing disks and knots')
    fix_disk_knots(input_disk, input_knots, output_disk, output_knots)
    print('Processing bulges')
    fix_bulge(input_bulge, output_bulge)
    print('Gzipping....')
    os.system("gzip -f %s %s %s"%(output_disk,output_bulge,output_knots))
    print('Done.')
