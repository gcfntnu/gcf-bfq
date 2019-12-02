import configparser
import os

def getConfig() :
    config = configparser.ConfigParser()
    config.read_file(open("/config/bcl2fastq.ini","r"))
    if("Paths" in config.sections()) :
        return config
    return None
