#!/usr/bin/env python
"""
Script to generate Run3.0i instance catalogs using scripts from the
SLSprinkler package.
"""
import os
import sys
import json
import argparse
import subprocess
import multiprocessing

class GenerateInstcat:
    """
    Class to run command line instance catalog generation scripts.
    """
    def __init__(self, config_file):
        self.config_file = config_file
        with open(config_file) as _:
            self.config = json.load(_)
        self.instcat_dir = self.config['instcat_dir']
        self.fov = self.config['fov']
        slsprinkler_dir = self.config['slsprinkler_dir']
        self.script_dir = os.path.join(slsprinkler_dir, 'scripts', 'dc2')
        if not os.path.islink('data'):
            os.symlink(os.path.join(slsprinkler_dir, 'data'), 'data')
        self.opsim_db_file = self.config['db']

    @staticmethod
    def _subprocess_check_call(command, dry_run):
        print(command)
        print()
        sys.stdout.flush()
        if not dry_run:
            subprocess.check_call(command, shell=True)

    def __call__(self, visit, dry_run=False):
        commands = []
        # Generate instcat for un-sprinkled objects.
        commands.append('time generateInstCat.py '
                        f'--config_file {self.config_file} '
                        f'--fov {self.fov} --ids {visit} --agn_threads 1 '
                        f'--out_dir {self.instcat_dir}')

        # Instcats for sprinkled objects:
        visit_dir = f'{visit:08d}'
        phosim_cat_file = os.path.join(self.instcat_dir, visit_dir,
                                       f'phosim_cat_{visit}.txt')

        # Strongly lensed AGNs
        agn_truth_cat = os.path.join(self.config['truth_table_dir'],
                                     self.config['agn_truth_cat_basename'])
        file_out = os.path.join(self.instcat_dir, visit_dir,
                                f'lensed_agn_{visit}.txt')
        commands.append(f'time python {self.script_dir}/create_agn_ic.py '
                        f'--obs_db {self.opsim_db_file} '
                        f'--obs_id {visit} '
                        f'--agn_truth_cat {agn_truth_cat} '
                        f'--file_out {file_out}')
        commands.append(f'gzip -f {file_out}')
        commands.append(f'echo "includeobj {os.path.basename(file_out)}.gz" '
                        f'>> {phosim_cat_file}')
        commands.append('echo')

        # Strongly lensed SNe
        sne_truth_cat = os.path.join(self.config['truth_table_dir'],
                                     self.config['sne_truth_cat_basename'])
        output_dir = os.path.join(self.instcat_dir, visit_dir)
        cat_file_name = f'lensed_sne_{visit}.txt'
        file_out = os.path.join(output_dir, cat_file_name)
        commands.append(f'time python {self.script_dir}/create_sne_ic.py '
                        f'--obs_db {self.opsim_db_file} '
                        f'--obs_id {visit} '
                        f'--sne_truth_cat {sne_truth_cat} '
                        f'--output_dir {output_dir} '
                        '--sed_folder Dynamic '
                        f'--cat_file_name {cat_file_name}')
        commands.append(f'gzip -f {file_out}')
        commands.append(f'echo "includeobj {os.path.basename(file_out)}.gz" '
                        f'>> {phosim_cat_file}')
        commands.append('echo')

        # Lensed host galaxies
        host_truth_cat = os.path.join(self.config['truth_table_dir'],
                                      self.config['host_truth_cat_basename'])
        fits_stamp_dir = self.config['fits_stamp_dir']
        file_out = os.path.join(self.instcat_dir, visit_dir,
                                f'lensed_hosts_{visit}.txt')
        commands.append('time python '
                        f'{self.script_dir}/create_lensed_host_ic.py '
                        f'--obs_db {self.opsim_db_file} '
                        f'--obs_id {visit} '
                        f'--fov {self.fov} '
                        f'--host_truth_cat {host_truth_cat} '
                        f'--fits_stamp_dir {fits_stamp_dir} '
                        f'--file_out {file_out}')
        commands.append(f'gzip -f {file_out}')
        commands.append(f'echo "includeobj {os.path.basename(file_out)}.gz" '
                        f'>> {phosim_cat_file}')

        command = ';\n'.join(commands)
        self._subprocess_check_call(command, dry_run)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Script to generate Run3.0i instance catalogs.')
    parser.add_argument('config_file', type=str, help='config file.')
    parser.add_argument('visit', type=str, help='visit to process or '
                        'text file with a list of visits.')
    parser.add_argument('--dry_run', default=False, action='store_true',
                        help='flag to do a dry run.')
    parser.add_argument('--processes', type=int, default=1,
                        help='number of processes to use.')
    parser.add_argument('--start', type=int, default=0,
                        help='starting line in text file with visits.')
    parser.add_argument('--end', type=int, default=None,
                        help='ending line in text file with visits.')
    args = parser.parse_args()

    try:
        visits = [int(args.visit)]
    except ValueError:
        with open(args.visit) as fd:
            visits = [int(_.strip().split()[0]) for _ in fd
                      if not _.startswith('#')]
        end = len(visits) if args.end is None else args.end
        visits = visits[args.start: end]

    processes = min(args.processes, len(visits))
    gen_instcat = GenerateInstcat(args.config_file)
    if processes == 1:
        for visit in visits:
            gen_instcat(visit, args.dry_run)
    else:
        with multiprocessing.Pool(processes=processes) as pool:
            workers = []
            for visit in visits:
                workers.append(pool.apply_async(gen_instcat,
                                                (visit, args.dry_run)))
            pool.close()
            pool.join()

            _ = [worker.get() for worker in workers]
