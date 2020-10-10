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
    maxVersion = allowedVersions[len(allowedVersions)-1]
    minVersion = allowedVersions[0]
    version = importlib.import_module(o.package).__version__

    subversions = version.split(".")
    newerVersion = True
    for i in range(0,len(subversions)):
        if int(subversions[i]) >= int(maxVersion.split(".")[i]) and newerVersion==True:
            newerVersion = True
        else:
            newerVersion = False

    if version in allowedVersions:
        print('Package %s version %s is compatible.' % (o.package,version))
    elif newerVersion == True:
        print('Package %s version %s seems to be newer than the latest confirmed version number %s. It may or may not work. Good luck!' % (o.package,version,maxVersion))
    else:
        print('Warning: package %s version %s does not seem to be compatible. Please upgrade to version %s or higher.' % (o.package,version,maxVersion))
