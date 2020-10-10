#!/usr/bin/env python
## GJH
import importlib, yaml, argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--package',required=True)
    parser.add_argument('-y', '--yamlfile',required=True)
    o = parser.parse_args()

    ## yaml file handling and extraction
    with open(o.yamlfile, 'r') as stream:
        try:
            yamldata=yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    allowedVersions = yamldata[o.package].split(',')

    version = importlib.import_module(o.package).__version__

    if version in allowedVersions:
        print('Package %s version %s is compatible.' % (o.package,version))
    else:
        print('Warning: package %s version %s is not found in the list of confirmed compatible versions. The highest version number that was confirmed to work is %s.' % (o.package,version,allowedVersions[len(allowedVersions)-1]))
