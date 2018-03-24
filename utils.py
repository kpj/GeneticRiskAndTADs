import yaml


def load_config(fname='config.yaml'):
    with open(fname) as fd:
        return yaml.load(fd)
