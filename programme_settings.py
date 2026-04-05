#############################################################
# See https://docs.python.org/3/library/configparser.html   #
# Holds the main settings for the applications              #
#############################################################

import configparser

settings = configparser.ConfigParser(allow_no_value=True)

def read():
    global settings
    with open("settings.ini") as fh:
        settings.read_file(fh)


def write():
    global settings

    with open("settings.ini", 'w') as configfile:
        settings.write(configfile)