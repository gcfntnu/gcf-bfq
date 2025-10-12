import configparser


def getConfig():
    config = configparser.ConfigParser()
    config.read_file(open("/config/bcl2fastq.ini"))
    if "Paths" in config.sections():
        return config
    return None
