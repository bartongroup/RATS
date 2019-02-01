#!/homes/kfroussios/bin/python3

"""mylogs.py

Author: Kimon Froussios
Last revised: 27/10/2017

Library for custom logging.

XXXstring() functions return formatted strings for custom use by the caller.
log_XXX() functions format and record messages to log-files.

NOTE: XXXstring() return values are not intended to be passed as arguments to 
log_XXX(), as this will cause duplication and will mess up the formatting. Use 
XXXstring() to print info to the screen or within your output files. Use 
log_XXX() to record to log-files.
"""

import sys, datetime, os, subprocess


# Do everything in ONE write operation, to prevent race conditions when accessing the logfile. 


def tstamp():
    """Readable current date-time.

    Returns:
        (str)
    """
    return str(datetime.datetime.now())


def escapise(comnd):
    """Escape special characters in command.

    Any special characters would have arrived here by having been escaped in the command line.
    To log the command in a way it can be reused as-is, escapes must be re-added.

    Args:
        comnd[str]: List of strings to sanitise.
    Returns:
        [str]
    """
    for i,c in enumerate(comnd):
        # Escape existing single quotes, so I can surround strings with them without interference.
        if '\t' in c:
            # If tabs are needed, expansion of $'\t' requires the string to be unquoted, so characters must be escaped individually.
            comnd[i] = c.translate(str.maketrans({"\"": r"\"",
                                                  "\'": r"\'",
                                                  "$": r"\$",
                                                  "\t": r"$'\t'",
                                                  ">": r"\>",
                                                  "<": r"\<",
                                                  "|": r"\|",
                                                  "*": r"\*",
                                                  "#": r"\#",
                                                  "%": r"\%"}))
        # If no tabs are needed, single-quote if any special characters are present. (cleanest option)
        elif any((x in [' ' ,'"', '$', '@', ':', '^', '%', '*', '>', '<', '|', '#']) for x in c):
            c = c.replace("'", "\'")
            comnd[i] = "'" + c + "'"
        # Sometimes an empty string is an argument. That needs to be quoted for the log.
        elif c == '':
            comnd[i] = '\'\''
    return comnd


def paramstring(message=""):
    """Execution parameters log-string.
    
    I normally use this at the beginning of my output.
    
    Args:
        message(str): Optional message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n", " ")
    return "###INFO### " + tstamp() + "\t" + " ".join(
        sys.argv) + "\n###INFO### \t\t\t\t### CWD: " + os.getcwd() + "  ### PYTHON: " + sys.executable + (
               "\n###INFO### " if message else "") + message + "\n"


def donestring(message=""):
    """Done log-string.
    
    I normally use this at the end of tasks.
    
    Args:
        message(str): A custom message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n", " ")
    return "###INFO### " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + "Done " + message + ".\n"


def infostring(message=""):
    """Info log-string.
    
    I normally use this at the end of tasks.
    
    Args:
        message(str): A custom message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n", " ")
    return "###INFO### " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + message + "\n"


def warnstring(message=""):
    """Warning log-string.
    
    Args:
        message(str): The message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n", " ")
    return "#WARNING!# " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + message + "\n"


def errstring(message=""):
    """Error log-string.
    
    Args:
        message(str): The message to add.
    Returns:
        (str)
    """
    message.rstrip().replace("\n", " ")
    return "##!ERROR!# " + tstamp() + "\t" + os.path.basename(sys.argv[0]) + " - " + message + "\n"


def log_command(message="", logfile="./commands.log"):
    """Record timestamp, command-line call and optional message.
    
    This function obtains the command from sys.argv[].
    
    Args:
        message(str): Message to record in the log. (Default empty)
        logfile(str): File to write to (./commands.log).
    """
    cmnd = " ".join(sys.argv[:2]) + ' ' + " ".join(escapise(sys.argv[2:]))
    with open(logfile, 'a') as comlog:
        if message == "":
            comlog.write(tstamp() + "\t" + str(sys.executable) + " " + cmnd + "\n")
        else:
            comlog.write(tstamp() + "\t" + str(sys.executable) + " " + cmnd + "\n" + "                          \t" + message.rstrip().replace("\n", "\n          ") + "\n")


def log_message(message="", logfile="./messages.log"):
    """Record timestamp and message.
    
    Records timestamp and message to specified log-file.
    
    Args:
        message(str): Message to record in the log. (Default empty)
        logfile(str): File to write to. (./messages.log)
    """
    with open(logfile, 'a') as comlog:
        comlog.write(tstamp() + "\t" + message.rstrip().replace("\n", "\n          ") + "\n")


######################
##### EXECUTABLE #####
######################


if __name__ == "__main__":

    if sys.argv[1] == "-e":
        # Log command and run it.
        c = " ".join(sys.argv[2:])
        log_command(message=c)  # Log both the full mylogs command (timestamped), and the command actually executed (as the message).
        subprocess.call(c, shell=True, stdout=sys.stdout)
    elif sys.argv[1] == "-m":
        # Log message.
        log_message(logfile=sys.argv[2], message=" ".join(escapise(sys.argv[3:])))

    sys.exit(0)

# EOF
