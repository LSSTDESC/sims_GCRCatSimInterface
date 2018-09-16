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

def main(in_instcat_disk, in_instcat_bulge, in_instcat_knots,
         out_instcat_disk, out_instcat_bulge, out_instcat_knots):

    # Use .fopen to read in the command and object lines from the
    # instance catalog.
    count = 0
    skipped = 0
    with fopen(in_instcat_disk, mode='rt') as input_disk,   \
         fopen(in_instcat_bulge, mode='rt') as input_bulge,  \
         fopen(in_instcat_knots, mode='rt') as input_knots,  \
         open(out_instcat_disk, 'w') as output_disk, \
         open(out_instcat_bulge, 'w') as output_bulge, \
         open(out_instcat_knots, 'w') as output_knots:

        # We first go through the knots catalog because some entries are missing
        # compared to the full bulge/disk catalog (faint knots have already been
        # removed at the instance catalog creation level)
        for line_knots in input_knots:

            # Extract the galaxy ID for that knots component
            tokens_knots = line_knots.strip().split()
            id_knots = int(tokens_knots[1]) >> 10

            found = False
            offensive = False
            # Loop through the disk and bulge catalogs
            for line_disk, line_bulge in zip(input_disk, input_bulge):
                tokens_disk = line_disk.strip().split()
                id_disk = int(tokens_disk[1]) >> 10

                tokens_bulge = line_bulge.strip().split()
                id_bulge = int(tokens_bulge[1]) >> 10

                # Checking that both components are for the same object
                if id_disk != id_bulge:
                    print("ERROR: Mismatch between disk and bulge catalogs", id_disk, id_bulge)
                    exit(-1)

                # Extracting the extinction parameters:
                if tokens_disk[17].lower() != 'none':
                    disk_internal_av = float(tokens_disk[18])
                    disk_internal_rv = float(tokens_disk[19])
                else:
                    disk_internal_av = 0
                    disk_internal_rv = 0

                if tokens_bulge[17].lower() != 'none':
                    bulge_internal_av = float(tokens_bulge[18])
                    bulge_internal_rv = float(tokens_bulge[19])
                else:
                    bulge_internal_av = 0
                    bulge_internal_rv = 0


                # Applying check on extinction parameters
                if ((disk_internal_av < 0) or (bulge_internal_av < 0) or
                    (disk_internal_rv > 1) or (bulge_internal_rv > 1)):

                    # Detected offensive object that shall be unceremoniously dropped
                    print('Dropping offensive galaxy %d with disk/bulge av/rv %f/%f %f/%f'%(id_disk,
                            disk_internal_av,disk_internal_rv,bulge_internal_av,bulge_internal_rv))
                    offensive = True
                    skipped +=1
                    # We are not printing these entries in the output instance catalog,
                    # and will also skip that knots entry
                    break

                if id_disk == id_knots:
                    found=True
                    break
                else:
                    output_disk.write(line_disk.strip()+'\n')
                    output_bulge.write(line_bulge.strip()+'\n')

            if offensive:
                # Ignoring that knots entry because it was found to be offensive
                continue

            if not found:
                print("ERROR: object ids do not match between input catalogs")
                exit(-1)

            # Get total flux
            magnorm_disk = np.float(tokens_disk[4])
            magnorm_knots = np.float(tokens_knots[4])
            total_flux = 10.**(-magnorm_disk/2.5) + 10.**(-magnorm_knots/2.5)
            knots_flux_ratio = 10.**(-magnorm_knots/2.5) / total_flux

            # Apply flux cap for large galaxies
            size = np.float(tokens_disk[13])
            if size > 2.5:
                knots_flux_ratio = np.clip(knots_flux_ratio, 0, 0.3)
                count+=1
                print("Capping knots flux for object %d, with magnorm: %f and size %f"%(id_knots,magnorm_disk,size))
            elif size > 1.:
                knots_flux_ratio = np.clip(knots_flux_ratio, 0, 0.5)
                count+=1
                print("Capping knots flux for object %d, with magnorm: %f and size %f"%(id_knots,magnorm_disk,size))

            magnorm_disk = -2.5*np.log10((1-knots_flux_ratio)*total_flux)
            magnorm_knots = -2.5*np.log10(knots_flux_ratio*total_flux)

            # Update the entry
            tokens_disk[4] = ("%.7f"%magnorm_disk).rstrip('0')
            tokens_knots[4] = ("%.7f"%magnorm_knots).rstrip('0')
            line_disk = ' '.join(tokens_disk)
            line_knots = ' '.join(tokens_knots)

            # Write
            output_disk.write(line_disk.strip()+'\n')
            output_knots.write(line_knots.strip()+'\n')
            output_bulge.write(line_bulge.strip()+'\n')

    print("Fixed %d galaxies"%count)
    print("Skipped %d galaxies"%skipped)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Knots cancelling script')
    parser.add_argument('input_disk', type=str)
    parser.add_argument('input_bulge', type=str)
    parser.add_argument('input_knots', type=str)
    parser.add_argument('output_disk', type=str)
    parser.add_argument('output_bulge', type=str)
    parser.add_argument('output_knots', type=str)
    args = parser.parse_args()
    main(args.input_disk, args.input_bulge, args.input_knots,
         args.output_disk,args.output_bulge, args.output_knots)
