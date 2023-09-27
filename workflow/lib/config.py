import os
import re
import sys
import csv

import collections.abc
from pathlib import Path
from Bio import SeqIO


def load_sample_list(samplelist_fp, paired_end, root_proj = ''):
    """
    Build a list of samples from a sample list file.
    :param samplelist_fp: a Path to a whitespace-delimited samplelist file,
       where the first entry is the sample name and the rest is ignored.
    :returns: A dictionary of samples with sample name and associated file(s)
    """
    Samples = {}
    with open(str(samplelist_fp)) as f:
        reader = csv.DictReader(f, fieldnames=['sample','1','2'])
        for row in reader:
            sample = row['sample']
            try:
                r1 = _verify_path(row['1'].rstrip())
            except ValueError:
                raise ValueError("Associated file for {} not found.".format(sample))
            r2 = ''
            if paired_end:
                try:
                    r2 = _verify_path(row['2'].rstrip())
                except ValueError:
                    raise ValueError("Paired-end files specified, but mate pair for '{}' is missing or does not exist.".format(sample))
            Samples[sample] = {'1': r1, '2': r2}
    return Samples


def read_seq_ids(fasta_fp):
    """
    Return the sequence identifiers for a given fasta filename.
    """
    return [record.id for record in SeqIO.parse(str(fasta_fp), "fasta")]


def makepath(path):
    return Path(path).expanduser()


def _verify_path(fp):
    if not fp:
        raise ValueError("Missing filename")
    path = Path(fp)
    if not path.is_file():
        raise ValueError("File not found")
    return str(path.resolve())


def verify(path):
    path = Path(path)
    if path.exists():
        return path.resolve()
    else:
        raise ValueError("Path %s does not exist" % path)


def validate_paths(cfg, root):
    """Process paths in config file subsection.

    For each key ending in _fp, the value is:
    - converted to a pathlib.Path
    - ensured to be an absolute path, by appending `root` if needed
    - ensured to exist, if it is not the value from `output_fp`
    - expanded home directory ~

    :param cfg: a config file subsection
    :returns: an updated copy of cfg
    """
    new_cfg = dict()
    for k, v in cfg.items():
        if k.endswith('_fp'):
            v = makepath(v)
            if not v.is_absolute():
                v = root/v
            if k != 'output_fp':
                try: v = verify(v)
                except ValueError:
                    raise ValueError(
                        "For key '%s': path '%s' does not exist" % (k,v))
        new_cfg[k] = v
    return new_cfg


def validate_config(cfg):
    """Resolve root in config file, then validate paths."""
    if 'root' in cfg['all']:
        root = verify(cfg['all']['root'])
    else:
        root = Path.cwd()
    # Iteratively check paths for each subsection
    new_cfg = dict()
    for section, values in cfg.items():
        if isinstance(values, dict):
            new_cfg[section] = validate_paths(values, root)
        else:
            new_cfg[section] = values
    new_cfg['all']['root'] = root
    return new_cfg


def output_subdir(cfg, section):
    return cfg['all']['output_fp']/cfg[section]['suffix']


def process_databases(db_dict):
    """Process the list of databases.

    Expands the nucleotide and protein databases specified
    """
    dbs = {'nucl':{}, 'prot':{}}
    root = verify(makepath(db_dict['root_fp']))
    nucl = db_dict.get('nucleotide')
    prot = db_dict.get('protein')
    if nucl:
        dbs['nucl'] = {db: str(root/path) for db, path in nucl.items()}
    if prot:
        dbs['prot'] = {db: str(root/path) for db, path in prot.items()}
    return dbs


def _update_dict(target, new):
    for k, v in new.items():
        if isinstance(v, collections.abc.Mapping):
            # We could use .get() here but ruamel.yaml's weird Mapping
            # subclass outputs errors to stdout if the key doesn't exist
            if k in target:
                target[k] = _update_dict(target[k], v)
            else:
                target[k] = _update_dict({}, v)
        else:
            target[k] = v
    return target


def _update_dict_strict(target, new):
    for k, v in new.items():
        if isinstance(v, collections.abc.Mapping) and k in target.keys():
            target[k] = _update_dict_strict(target.get(k, {}), v)
        elif k in target.keys():
            target[k] = v
        else:
            sys.stderr.write("Key '%s' not found in target, skipping\n" % k)
            continue
    return target


def update(config_str, new, strict=False):
    config = ruamel.yaml.round_trip_load(config_str)
    if strict:
        config = _update_dict_strict(config, new)
    else:
        config = _update_dict(config, new)
    return config
