#!/homes/kfroussios/bin/python3

"""fileutilities.py

Author: Kimon Froussios
Compatibility tested: python 3.5.2
Last reviewed: 13/01/2017

This module is a solution for Frequently Performed Generic Tasks that involve 
multiple files:
* repeating a command for a range of files (not currently parallelized),
* accessing and restructuring (multiple) delimited files.
* miscellaneous stuff. Some of it is auxiliary to the primary functions, some 
  is a legacy of this module's evolution of concept. Some of it could arguably 
  be in different libraries, but I see not reason to inflate the dependencies. 
  I hate sprawling webs of dependencies on non-standard packages more than I 
  dislike some out-of-place functions.

The primary purpose of this entire module is to encapsulate frequent generic 
tasks and related boilerplate code so as to be able to do them from the command 
line in a single step, instead of implementing specialized code for each 
variation of the task. The module provides a library of flexible functions as 
well as a main() implementing the primary use scenarios.

Execute with -h in a shell to obtain syntax and help.
"""

# This module consists of:
# -    a class for handling lists of files,
# -    a library of functions that perform generic tasks on multiple files, and
# -    a main that provides access to most of the above functionality

# NOTE about DataFrame indexes and headers:
# Although dataframes support row and column labels, these make content manipulations 
# in the context of this module harder. Instead, any labels present in the text 
# input are treated as plain rows or columns. These may be optionally dropped or 
# preserved, but the dataframe in-built labels for columns and rows are reserved 
# solely for custom use. When appropriate these custom labels will be included in 
# the output.

####### TODO SUGGESTIONS #######
#
#
################################

######## UPDATES ###############
#
# 2016-04-13 : Range support added to --cols. Closed-end ranges only.
# 2016-04-13 : Get columns by name support added to --cols.
# 2016-04-13 : Multiple ranges support added for --loop R. NEW SYNTAX!
# 2016-04-13 : Comma separated value support added to TARGET and --cols ONLY! Old syntax also retained.
#              Comma lists are not practical for other flags as they accept regex strings. 
# 2016-04-13 : Added option to retain fixed number of metadata lines at the top of the files
#              --cols, --rndcols and --appnd.
#
################################

import os, sys, string, re, subprocess, random, argparse
import pandas as pd
from builtins import list
from collections import Counter

import mylogs as ml


#####   F U N C T I O N S   #####


# http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
def natural_sorted(l):
    """Sort list of numbers/strings in human-friendly order.
    
    Args:
        l(list): A list of strings.
    Returns:
        list
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


def expand_fpaths(flist):
        """Fully expand and absolute-ify the paths of listed files.
        
        Does not verify path validity. All paths are expanded.
        
        Args:
            flist[str]: A list/FilesList of files.
        Returns:
            [str]: List of expanded paths.
        """
        return [os.path.abspath(os.path.expanduser(str(f))) for f in flist]


def are_empty(flist, invalids=True):   
    """List the files that are empty.
    
    Args:
        flist: A list/FilesList of files to probe.
        invalid(bool): How should invalid paths be reported? (Default True)
                    There is no right or wrong way to classify invalid
                    paths, it is entirely dependent on the context in which
                    the files would be used downstream. 
    Returns:
        FilesList: List of empty files.
    """
    try:
        flist.aliases[0]
    except AttributeError:
        # If it's a plain list, upgrade it to a FilesList.
        flist = FilesList(flist)
    result = FilesList()
    for i, (myfile, myalias) in flist.enum():
        try:
            if os.path.getsize(myfile) == 0:
                result.append(myfile, myalias)
        except OSError:
            if invalids:
                result.append(myfile, myalias)
    return result


def dont_exist(flist): 
    """List the files that don't exist.
    
    Args:
        flist: A list/FilesList of files to probe.
    Returns:
        FilesList: List of missing files.
    """
    try:
        flist.aliases[0]
    except AttributeError:
        # If it's a plain list, upgrade it to a FilesList.
        flist = FilesList(flist)
    result = FilesList()
    for i, (myfile, myalias) in flist.enum():
        if not os.path.exists(myfile):
            result.append(myfile, myalias)
    return result
    

def arent_readable(flist): 
    """List inaccessible files.
    
    Invalid paths are also reported as inaccessible.
    
    Args:
        flist: A list/FilesList of files to probe.
    Returns:
        FilesList: List of inaccessible files.
    """
    try:
        flist.aliases[0]
    except AttributeError:
        # If it's a plain list, upgrade it to a FilesList.
        flist = FilesList(flist)
    result = FilesList()
    for f, (myfile, myalias) in flist.enum():
        try:
            if not os.access(myfile, os.R_OK):
                result.append(myfile, myalias)
        except OSError:
            result.append(myfile, myalias)
    return result


def arent_text(flist, invalids=True):
    """List the file that are (probably) text.
    
    This is achieved by probing the contents and using heuristics to
    determine the type of those contents.
    
    Args:
        flist: A list/FilesList of files to probe.
        invalids(bool): Should invalid paths be reported as not-text?
                    (Default True)
    Returns:
        FilesList: List of files that are probably not plain text.
    """
    try:
        flist.aliases[0]
    except AttributeError:
        # If it's a plain list, upgrade it to a FilesList.
        flist = FilesList(flist)
    result = FilesList()
    for i, (myfile, myalias) in flist.enum():
        try:
            with open(myfile) as f:
                if not istext(f.read(1024)):
                    result.append(myfile, myalias)
        except IOError:
            if invalids:
                result.append(myfile, myalias)
    return result

# Helper function, string check.
def istext(s):
    """Check if a string is (probably) text.
    
    Use heuristic based on characters contained in the string, adapted from:
    http://code.activestate.com/recipes/173220-test-if-a-file-or-string-is-text-or-binary/
    
    Args:
       s(str): A string to test.
    Returns:
        bool
    """
    # Copy-pasted. No idea what the code means.
    text_characters = "".join(list(map(chr, list(range(32, 127)))) + list("\n\r\t\b"))
    _null_trans = string.maketrans("", "")
    if "\0" in s:
       return False
    if not s:  # Empty files/strings are considered text
        return True
    # Get the non-text characters (maps a character to itself then
    # use the 'remove' option to get rid of the text characters.)
    t = s.translate(_null_trans, text_characters)
    # If more than 30% non-text characters, then
    # this is considered a binary file
    if float(len(t))/float(len(s)) > 0.30:
        return False
    return True


def slink(flist, aliases=None, dir="./", autoext=True):
    """Create symbolic links for multiple files.
    
    Create a link for each of the listed paths into the specified directory,
    using the specified aliases. Items in the lists will be matched one for 
    one.
    If the aliases argument is omitted, the names for the links will be drawn
    from the aliases attribute of the paths list, if it is a FilesList object.
    If no aliases exist in either form, the files will be linked in the current
    or specified directory, using names their current basename.
    
    If linking to files of the same name located in different directories, a
    number will be automatically suffixed to the basename.
    
    Args:
        flist[str]: A list/FilesList of paths to link to.
        aliases[str]: A list of respective names for the created links. If 
                    omitted, the alias attribute of the flist argument will be 
                    used, and failing that, the existing basenames will be used.
        dir(str): The path to the directory in which the links should be 
                    placed. (Default "./")
        autoext(bool): Add the file extensions to the created links, if the
                    links are created from aliases that lack them.
                    (Default True)
    """
    if not aliases:
        # No link names provided. Try to find them elsewhere or create them.
        try:
            # flist is a FilesList and has the aliases attribute.
            aliases = flist.aliases
        except AttributeError:
            # flist is a plain list, so compute the link name from the file name.
            aliases = [os.path.basename(p) for p in flist]
    # Check for duplicate aliases and amend them. 
    # This applies mainly to link names automatically created from filenames, as 
    # the same file name can exist in different directories.
    if len(set(aliases)) < len(flist):
        aliases = autonumerate(aliases)
    # Add extensions where necessary, if desired.
    if autoext:
        for i in range(0, len(flist)):
            (b, p) = os.path.splitext(flist[i])
            c = p
            # If it's a .gz, include the next nested extension as well.
            if p == ".gz":
                p = os.path.splitext(b)[1] + p
            # Don't duplicate the extension if the alias already has it.
            a = os.path.splitext(aliases[i])[1]
            if c != a:
                aliases[i] = aliases[i] + p
    # Link.
    for i, mypath in enumerate(flist):
        os.symlink(mypath, os.path.join(dir, aliases[i]))


# Helper function.
def autonumerate(things):
    """Detect duplicate entries in a string list and suffix them.
    
    Suffixes are in _N format where N a natural number >=2. Existing suffixes
    in that format will also be detected and incremented.
    
    Args:
        things[str]: A list of strings.
    Returns:
        [str]: A corrected list of strings.
    """
    c = Counter(things);
    # Because I use decrement, reversing the list ensures first instance gets smallest number.
    things.reverse()
    for i, t in enumerate(things):
        n = c[t]
        if n > 1:  # The first occurrence is not suffixed. 
            newname = t +'_' + str(n)
            while newname in things:  # Check for already present suffixes
                n += 1
                newname = t +'_' + str(n)
            things[i] = newname
            c[t] -= 1
    things.reverse()
    return things


def make_names(items, parameters):
    """Automatically create file names based on parameters.
    
    If automatic names happen to turn out identical with one another, unique
    numbers are appended to differentiate them. Check documentation for
    autonumerate().
    
    Args:
        items[str]: A list of strings/filenames/paths to use as the basis for
                    the output names.
        parameters(str,str,str): The first element is the output directory, 
                    the second is a common prefix to add to the names, 
                    the third is a common sufix to add to the names.
                    Like so: <out[0]>/<out[1]>item<out[2] .
                    If any of the 3 values in None, no outnames will be made.
                    Use current directory and empty strings as necessary.
    Returns:
        [str]: A list of file paths.
    """
    outfiles = []
    if None not in parameters:
        for i in items:
            outfiles.append(os.path.join(os.path.abspath(os.path.expanduser(parameters[0])), 
                                         parameters[1] + os.path.splitext(os.path.basename(str(i)))[0] + parameters[2]) )
        autonumerate(outfiles)
    return outfiles


def do_foreach(flist, comm, comments=False, progress=True, out=(None,None,None), log=False):
    """Execute an arbitrary command for each of the listed files.
    
    Enables executing a shell command over a range of items, by inserting the
    item values into the command as directed by place-holder substrings.
    Although the above is how it is meant to be used, the values in the
    FilesList could be overridden to be any arbitrary string, in which case
    only two of the placeholders will work correctly (full and alias). The
    others are computed from these and may not be sensible when the items are
    not file paths.
    
    This is the only function with comments or progress attributes, because 
    invoked commands print their own output directly, so any informative messages 
    controlled by this library need to also be inserted in real time.
    
    Args:
        flist[]: A FilesList. If a plain list is given, aliases will be 
                    automatically computed.
        comm[str]: The components of an arbitrary command, with place-holders
                    for the item of each iteration. Placeholder strings are:
                    '***full***' : absolute path of file.
                    '***path***' : absolute path to the file.
                    '***file***' : filename without path.
                    '***core***' : filename without the last extension.
                    '***alias***': alias for file, by default same as core.
        comments(bool): Print commented call details to STDOUT. (Default False)
        progress(bool): Show start and completion of iterations on STDERR. 
                    (Default True)
        out(str,str,str): The first element is the output directory, the second
                    is a common prefix to add to the names, the third is a 
                    common suffix to add to the names. Check documentation for
                    make_names().
        log(bool): Log to /commands.log each individual call.
    """
    outstream= sys.stdout
    # Test if it is a FilesList or plain list. Upgrade it if it's plain. 
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Create output files. [] if out contains None.
    outfiles = make_names(flist, out)
    for i, (myfile, myalias) in flist.enum():
        # Substitute place-holders.
        command = []
        for c in comm:
            (mypath, mybase) = os.path.split(str(myfile))
            c = c.replace("***full***", str(myfile))
            c = c.replace("***path***", mypath)
            c = c.replace("***file***", mybase)
            c = c.replace("***core***", os.path.splitext(mybase)[0])
            c = c.replace("***alias***", str(myalias))
            command.append(c)
        # Redirect output.
        if outfiles:
            outstream = open(outfiles[i], 'w')
        # Do.
        try:
            see = " ".join(command)
            if log:
                ml.log_message(message=see, logfile="./commands.log")
            if comments and out == (None,None):
                outstream.write(ml.infostring("CWD: "+ os.getcwd() +"\tDO: "+ see))
            if progress:
                sys.stderr.write(ml.infostring("DO: "+ see))
        except IOError:
            pass
        subprocess.call(command, stdout=outstream, shell=False)
        # Optionally identify iteration.
        try:
            if comments and out == (None,None):
                outstream.write(ml.infostring("Finished: "+ str(myalias) +"\n"))
            if progress:
                sys.stderr.write(ml.infostring("Finished: "+ str(myalias) +"\n"))
        except IOError:
            pass
        finally:
            if outfiles:
                outstream.close()
        

def tails(flist, linesToGet=5):
    """Get last lines of text files.
     
    Adapted from:
        https://gist.github.com/volker48/3437288
     
    Files that cause an error are still represented in the return value with an 
    empty sublist, like empty files, in order to maintain the correspondence
    between input and output lists.
         
    Args:
        flist: A list/FilesList of text files.
        lines_to_get: Number of lines to get from the end.
    Returns:
        [[str]]: A list of lists. Each sublist represents the lines returned by
                    one file.
    Raises:
        IOError: With a message about ALL the inaccessible files.
    """
    result = []
    errors = []
    for myfile in flist:
        if linesToGet < 1:
            result.append([])
            continue
        try:
            with open(myfile, 'r') as f:
                # Place the cursor at the last character.
                f.seek(-1, os.SEEK_END)
                position = f.tell()
                linesSeen = 0
                # I want the beginning of a line, by proxy of the previous line's newline end.
                # So if the last character of the file is a newline, ignore it and go one step back.
                if f.read(1) == "\n":
                    position -= 1
                    f.seek(position)
                # Crawl backwards one character at a time.
                while linesSeen < linesToGet and position > 0:
                    # Try previous character.
                    position -= 1
                    f.seek(position)
                    c = f.read(1)
                    if c == "\n":
                        linesSeen += 1
                if position == 0:
                    # The cursor will be right after a newline, unless I hit the beginning
                    # of the file, in which case the loop leaves it on the second character.
                    # So bring it manually to the start of the line and file. 
                    f.seek(0)
                # Now that the cursor is a the right place, read in the lines.
                chunk = []
                for line in f:
                    chunk.append(line)
                result.append(chunk)
        except IOError as e:
            result.append([])
            errors.append(e.message)
    if errors != []:
        raise IOError(" --- ".join(errors))
    return result
 
 
def heads(flist, linesToGet=5):
    """Get first lines of text files.
     
    Files that cause an error are still represented in the return value with an 
    empty sublist, like empty files, in order to maintain the correspondence
    between input and output lists.
     
    Args:
        flist: A list/FilesList of text files.
        lines_to_get: Number of lines to get from the top.
    Returns:
        [[str]]: A list of lists. Each sublist represents the lines returned by
                    one file.
    Raises:
        IOError: With a message about ALL the inaccessible files. 
    """
    result = []
    errors = []
    for myfile in flist:
        if linesToGet < 1:
            result.append([])
            continue
        try:
            with open(myfile, 'r') as f:
                chunk = []
                for i in range(0, linesToGet):
                    line = f.readline()
                    if line != "":
                        # Otherwise it appends empty strings when it hits the end of file.
                        chunk.append(line)
                result.append(chunk)
        except IOError as e:
            result.append([])
            errors.append(e.message)
    if errors != []:
        raise IOError(" --- ".join(errors))
    return result        


def swap_strFiles(flist, insep=[","], outsep="\t"):
    """Replace the column separator with a different one.
    
    Supports multiple different delimiters in the input, to support one-step
    uniformity when the input files have different delimiters, but ALL input 
    will be split at ALL/ANY occurring delimiters. If the delimiter of one 
    file is present in a different use in an other file, the output may not 
    be what you want.
    Although made for converting delimited text, inseps and outsep could be any
    substring in a text, delimited or not.
    
    Args:
        flist: A list/FilesList of delimited text files.
        insep[str]: A list of regex strings. (Default [","])
        outsep(str): New column separator. (Default "\t")
    Returns:
        [str]: A list of strings with the changed delimiters. One string per 
                file. It is up to you to decide  what to do with the new 
                strings. The order of strings is the same as the input.
    """
    input = []
    if flist == []:
        # Read all data from STDIN at once. Input[] gets a single entry.
        input.append(sys.stdin.read())
    else:
        # Read in all the files entirely. Input[] gets as many entries as there are files.
        for myfile in flist:
            with open(myfile) as f:
                input.append(f.read())
    return swap_substr(input, insep, outsep)

# Helper function
def swap_substr(slist, insep=[","], outsep="\t"):
    """Replace all occurrences of insep with outsep.
    
    Insep may be a regex.
    
    Args:
        slist[str]: A list of strings.
        insep[str]: A list of regex strings. (Default [","])
        outsep(str): New substring. (Default "\t")
    Returns:
        [str]: A list of the edited strings. The order of the strings is the 
                    same as the input.
    """
    rx = re.compile("|".join(insep), re.MULTILINE)
    result = []
    for s in slist:
        # Replace all delimiters with the new one.
        result.append(rx.sub(outsep, s))
    return result

 
def prepare_df(df, myalias="", keyCol=None, keyhead="row_ID", header=False, cols=None, appendNum=True):
    """Prepare row names and column names.
    
    Assign column as row labels, rename columns based on their position and an
    arbitrary alias name for the dataframe, drop the first row.
    
    Args:
        df(pandas.DataFrame): A dataframe.
        myalias(str): The basename for the relabelling.
        header(bool): Remove first row (Default False).
        keyCol(int): Column to be used as row index. If None, no index will be
                    used. (Default None)
        keyhead(str): Label for the index.
        cols[int]: Custom index numbers for the columns (Default None). If None
                    then their current index positions are used. 
        appendNum(bool): Append the columns' positional indices to the alias
                    when making the new names (True).
    Returns:
        pandas.DataFrame
    """
    # Set row labels.
    if keyhead is None:
        keyhead = "row_ID"
    if keyCol is not None:
        # Add index without dropping it, so as not to affect column positions.
        df.set_index(df.columns.values.tolist()[keyCol], inplace=True, drop=False)
        df.index.name = str(keyhead)
    # Make custom column labels, based on alias and column position.
    if not cols:
        cols = list(range(0, df.shape[1]))
    labels = []
    if appendNum:
        labels = [str(myalias) +"_|"+ str(i) for i in cols]
    else:
        labels = [str(myalias) for i in cols]
    df.columns = labels
    # Remove header.
    if header:
        df.drop(df.index.values.tolist()[0], axis=0, inplace=True)
    return df
  

def count_columns(flist=[None], colSep=["\t"]):
    """Determine the number of fields in each file by inspecting the first row.
    
    Args:
        flist: A list of FilesList of files.
        colSep[str]: A list of characters used to separate columns.
    Returns:
        [int]: A list, in the same order as the given files.
    """
    tokenizer = re.compile("|".join(colSep))
    counts = []
    for file in flist:
        f = None
        if file is None:
            f = sys.stdin
            file = "<STDIN>"
        else:
            f = open(file)
        while True:
            line = f.readline()
            # Skip comments.
            if line[0] != "#":
                counts.append(len( tokenizer.split(line.rstrip()) ))
                break
            f.readline
        if f != sys.stdin:
            f.close()
    return counts


def get_valuesSet(flist=[None], axis='r', index=0, filter='a', colSep=["\t"]):
    """"List the set of different values in the column(s)/row(s).
    
    Args:
        flist: A list of FilesList of files.
        colSep[str]: A list of characters used to separate columns.
        index: Position index of the required column/row.
        axis(str): Data slice orientation - 'r' for row, 'c' for column.
        filter(str): non redundant set of: 'a' - all, 'u' - unique, 'r' - 
                    repeated values.
    Returns:
        [[]]: A list of lists. The inner lists represent the sets in order as 
                    requested.
    Raises:
        ValueError: Invalid axis or filter values.

    """
    tokenizer = "|".join(colSep)
    result = [] 
    if flist == []:
        # Use None as a flag to read from STDIN
        flist.append(None)        
    # Test if it is a FilesList or plain list. Upgrade it if it's plain. 
    # It will have at least one entry for sure by now, either way.
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Main part of this function.
    results = []
    for f, (myfile, myalias) in flist.enum():
        # Input.
        df = None
        instream = sys.stdin
        if myfile is not None:
            instream = open(myfile)
        df = pd.read_csv(instream, sep=tokenizer, header=None, index_col=None, comment="#", engine='python')            
        if instream != sys.stdin:
            instream.close()
        # Get value set.
        values = None
        if axis == 'r':
            values = df.iloc[int(index),:].tolist()
        elif axis == 'c':
            values = df.iloc[:,int(index)].tolist()
        else:
            raise ValueError("".join(["Unrecognized option: axis=", axis]))
        # Get count per value
        c = Counter(values);
        # Filter.
        if filter == 'a':
            results.append( set(values) )
        elif filter == 'u':
            results.append( set([v for v in values if c[v] == 1]) ) # set() is redundant but keeps output type consistent
        elif filter == 'r':
            results.append( set([v for v in values if c[v] > 1]) )
        else:
            raise ValueError("".join(["Unrecognized option: filter=", filter]))
    return results
    
    
def get_columns(flist=[None], cols=[0], colSep=["\t"], header=False, index=None, merge=True):
    """Obtain the specified columns.
        
    Comment lines starting with '#' are ignored.
    The data columns are assembled into a single DataFrame. 
    
    The returned columns will be labeled based on the name of the file they
    came from and their position in it. Existing labels are optionally 
    preserved as the top row or can be skipped entirely.
    
    If an index is specified, it will be used only for merging, and will NOT be
    included in the output columns, unless explicitly present in cols[].
    
    Args:
        flist: A list/FilesList of delimited plain text files.
        header(bool): Crop the header line (first non-comment line). (Default False)
        cols[int/str] : A list of positional indexes or names or ranges of the 
                    desired columns. (Default [0]).
        colSep[str]: List of characters used as field separators. 
                    (Default ["\t"]).
        merge(bool): Concatenate results from all files into a single 
                    dataframe. If False, a list of dataframes is returned
                    instead. (Default True).
        index(int): Column to be used as row index for merging. (Default None)
    Returns:
        [pandas.DataFrame]: List of DataFrames. If merge=True, only the 
                    first element will be populated.
    """
    tokenizer = "|".join(colSep)
    result = [] 
    if flist == []:
        # Use None as a flag to read from STDIN
        flist.append(None)        
    # Test if it is a FilesList or plain list. Upgrade it if it's plain. 
    # It will have at least one entry for sure by now, either way.
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Parse.
    keyhead = None
    for f, (myfile, myalias) in flist.enum():
        # I used to use the pandas parser, with my own parser used only as fallback
        # for problematic cases. As flexibility requirements increased, using the
        # pandas parser became too opaque and difficult to maintain, 
        # so now all cases are delegated to mine. 
        df = get_columns_manual(myfile, cols=cols, colSep=colSep, header=header, 
                                        alias=myalias, index=index)
        if not keyhead:
            keyhead = df.index.name
        result.append(df)
    # Merge.
    if merge:
        result = [pd.concat(result, axis=1, join='outer', ignore_index=False)]
        result[0].index.name = keyhead
    return result


# Helper function
def get_columns_manual(file=None, cols=[0], colSep=["\t"], header=False, index=None, alias=None):
    """Get specified columns from a file where rows have varying numbers of fields.
    
    Some tables contain a fixed set of fields followed by optional fields. In
    these rare cases, traditional parsers fail due to inconsistent number of 
    fields. This function provides a work-around for that.
    
    It is entirely the user's responsibility to ensure that the inconsistent 
    row lengths are not a symptom of table corruption/malformation and that it 
    is safe and reliable to extract the desired columns. If a row is shorter 
    than expected, it is padded with the value "IDXERROR". If this value shows 
    up in your result and you are not explicitly expecting it, you should stop
    and seriously examine your input table.
    
    Args:
        file(str): A delimited plain text file.
        header(bool): If True, the first non-comment line will not be in
                    the data. (Default False)
        cols[int]: A list of positional indexes of the desired columns. 
                    (Default [0]).
        colSep[str]: List of regex strings for field separators. 
                    (Default ["\t"]).
        index(int): Position of column to be used as row index. (Default None)
        alias(str): An alias for the file. Used for naming the columns.
    Returns:
        pandas.DataFrame: DataFrame with the columns, labeled by original 
                    column number, ordered as specified.
    """
    tokenizer = re.compile("|".join(colSep))
    # Input source.
    f = None
    if file is None:
        f = sys.stdin
        file = "STDIN"
    else:
        f = open(file)
    if alias is None:
        alias = FilesList.autoalias(file)
    # Import data.
    keyhead = None
    values = []
    labels = []
    for l, line in enumerate(f):
        if line[0] == '#' or line == "\n":
            # Skip comments and empty lines.
            continue
        else:
            # Get the fields.
            fields = tokenizer.split(line.rstrip("\n"))
            # Column labels from the first non-comment non-empty row, 
            # regardless of whether they really are labels or not.
            if not labels:
                labels = fields
            # Find out name of row index.
            if (not keyhead) and header and (index is not None):
                keyhead = str(fields[index])
            # Get columns.
            selection = []
            expandedcols = []
            for c in cols:
                v = str(c).split(":")
                if len(v) == 1:
                    try:
                        expandedcols.append(int(v[0]))
                    except ValueError:
                        expandedcols.append(labels.index(v[0]))
                else:
                    try:
                        expandedcols.extend(list(range(int(v[0]), int(v[1]) + 1)))
                    except TypeError:
                        expandedcols.extend(list(range(labels.index(v[0]), labels.index(v[1]) + 1)))
            for i in expandedcols:
                try:
                    selection.append(fields[i])
                except IndexError:
                    # Silently adding fields is too dangerous, so a flag value is needed. 
                    # Values like None or NA can sometimes be legitimate values for fields.
                    selection.append("IDXERROR")
            # Add the key at the end, where they won't interfere with column numbers.
            if index is not None:
                selection.append(fields[index])
            values.append(selection)
    if f != sys.stdin:
        f.close()
    # Adjust index of row keys to reflect the fact I stuck them at the end.
    if index is not None:
        index = len(values[0])-1
        expandedcols.append("my_garbage_label_row_key")
    # Package data nicely.
    df = pd.DataFrame(data=values)
    df = prepare_df(df, myalias=alias, keyCol=index, header=header, cols=expandedcols, 
                    keyhead=keyhead, appendNum=True if len(expandedcols)>1 else False)
    if index is not None:
        df.drop(alias+"_|my_garbage_label_row_key", 1, inplace=True)
    return df


def get_random_columns(flist, colSep=["\t"], k=1, header=False, index=None, merge=True):
    """ Get k random columns from each file.
    
    The returned columns will be labeled based on the name of the file they
    came from and their position in it. Existing labels are optionally 
    preserved as the top row or can be skipped entirely.
    
    If an index is specified, it will be used for merging (if applicable) and 
    will be included as a column in each output file.
    
    Args:
        flist: A list or FilesList of files.
        k(int): How many columns to get.
        colSep[str]: A list of characters used as field separators. 
                    (Default ["\t"])
        header(bool): Strip column headers. (Default False)
        index(int): Column to use as row index for merging. (Default None)
        merge(bool): Concatenate results from all files into a single 
                    dataframe. If False, a list of dataframes is returned
                    instead. (Default True).
    Returns:
        [pandas.DataFrame]: List of DataFrames. If merge=True, only the 
                        first element will be populated.
    """
    tokenizer = "|".join(colSep)
    # The files may have different number of columns
    fieldNums = count_columns(flist, colSep)
    result = []
    if flist == []:
        # Use None as a flag to read from STDIN
        flist.append(None)        
    keyhead = None
    # Test if it is a FilesList or plain list. Upgrade it if it's plain. 
    # get_columns() does this too, but as I call it per item in flist, I *must*
    # preserve any alias that is potentially already specified.
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Get columns.
    for f, (myfile, myalias) in flist.enum():
        cols = []
        if index is not None:
            # Generate random choice of columns.
            # range() is right-open. 
            cols.extend(random.sample(list(range(0, fieldNums[f]-1)), k))  
        else:
            cols = random.sample(list(range(0,fieldNums[f])), k)
        # Would normally delegate the actual getting to get_columns() but there 
        # are too many little differences to accommodate that complicate the code
        # to the point of defeating any benefits from code re-use. 
        df = pd.read_csv(myfile, sep=tokenizer, header=None, index_col=None, comment="#", engine='python')
        if (not keyhead) and header and (index is not None):    
            keyhead = str(df.iloc[0,index])
        # Adjust row and column labels.
        df = prepare_df(df, myalias=myalias, keyCol=index, header=header, keyhead=keyhead,
                        appendNum=True if k>1 else False)
        # Slice the part I need.
        df = df.iloc[:,cols] 
        result.append(df)
    # Merge. 
    if merge:
        result = [pd.concat(result, axis=1, join='outer', ignore_index=False)]
        result[0].index.name = keyhead
    return result
    

def append_columns(flist, colSep=["\t"], header=False, index=None):
    """Append all columns from the files, as they are.
    
    The order of the columns in the output follows the order of the columns in
    the files and the order of the files. The files don't need to have the same
    number of columns. It is your responsibility to ensure sensible consistency
    of number and order of rows across files, otherwise the combined data may
    be nonsensical.
    
    This function also supports key-aware appending, using outer-join, when a
    row index is specified.
    
    Args:
        flist: A list/FilesList of files to combine.
        colSep[str]: A list of characters used as field delimiters. 
                    (Default ["\t"])
        header(bool): First non-comment line as column labels. (Default False)
        index(int): Column to use as row index (same in all files). 
                    (Default None)
    Returns:
        pandas.Dataframe
    """
    try:
        flist.aliases[0]
    except AttributeError:
        flist = FilesList(flist)
    # Determine how many columns each file has.
    numofcols = count_columns(flist, colSep=colSep)
    # Delegate fetching all the columns.
    data = []
    keyhead = None
    for f, (myfile, myalias) in flist.enum():
        # List the columns and remove the index one from among them.
        cols = [i for i in range(0,numofcols[f]) if i != index]
        df =get_columns(FilesList(files=[myfile], aliases=[myalias]), cols=cols, 
                     colSep=colSep, header=header, merge=False, index=index)[0]
        data.append( df )
    # Merge. Row indexes will have been assigned by get_columns(), if applicable.
    keyhead = data[0].index.name
    result = pd.concat(data, axis=1, join="outer", ignore_index=False)
    result.index.name = keyhead
    return result


def get_crosspoints(flist, cols=[0], rows=[0], colSep=["\t"], header=False, index=None, merge=True):
    """ Get the values at selected rows and columns.
    
    The values at the intersections of the selected rows and columns are extracted.
    
    Args:
        flist: A [str] list or fileutilities.FilesList of delimited text files.
        colSep[str]: List of column separators.
        cols[int]: List of columns.
        rows[int]: List of rows.
        header(bool): Whether there is a header line (False).
        index(int): Which column has the row labels (None).
        merge(bool): Merge results into single table (True).
    Returns:
        [pandas.DataFrame]: 
    """
    results = get_columns(flist, cols=cols, colSep=colSep, header=header, merge=merge, index=index)
    for i in range(0, len(results)):
        results[i] = results[i].iloc[rows,:]
    return results
    
    
#####   C L A S S E S   #####



class FilesList(list):
    """A container for a list of files.

    An extension of the built-in list, specifically for files, providing a
    means to import multiple filenames either from text lists or from 
    directories. The purpose is to facilitate batch operations and sensible 
    output of their results.
    
    FilesList is generally backwards compatible with the built-in list and it
    should be possible for them to be used interchangeably in most cases. A 
    plain list can be cast as a FilesList, when necessary, allowing appointment
    of default alias values. A FilesList should always work as a plain list
    without additional actions (except for its string representation). When a 
    FilesList is accessed as a plain list, only the full paths will be 
    accessed. Certain equivalent methods are supplied for 
    
    Most basic operations inherited from list are supported. Appending has been
    overridden to keep paths and aliases in sync. Sorting, deleting and 
    inserting of items are currently not supported and will break the 
    correspondence between paths and aliases.
    
    Attributes defined here:
        aliases = [] : Practical aliases for the full file-paths.
    """
    def __init__(self, files=None, aliases=None, fromtuples=None):
        """Construct an instance of the FilesList.
        
        A FilesList can b created:
        - empty
        - from a list of files (with default aliases automatically assigned)
        - from a list of files and a list of aliases (in the same order)
        - from a list of (file, alias) tuples.
        
        Args:
            files[str]: A list of files. (Default None)
            aliases[str]: A list of aliases. (Default None)
            fromtuples[(str,str)]: A list of tuples (file, alias). (Default 
                        None) If this is specified together with flist and/or 
                        aliases, the data in fromtuples is used only.
        """
        # If data is a list of (file, alias) tuples, unpair tuples into lists.
        if fromtuples is not None:
            data = [list(t) for t in zip(*fromtuples)]
            # Any data passed to flist and aliases at method call is discarded.
            files = data[0]
            aliases = data[1]
        # Having aliases without the paths is rather useless.
        if aliases:
            if not files:
                raise ValueError("No files supplied for the aliases.")
        else:
            # Make None into empty.
            aliases = []
        # Assign default aliases to be same as files. Expand file paths.
        if files is not None:
            files = expand_fpaths(files)
            if not aliases:
                for f in files:
                    aliases.append(self.autoalias(f))
        else:
            # If still empty, it was an empty call to the constructor.
            files = []
        # Create the basic list.
        super(FilesList, self).__init__(files)
        # Add a plain list attribute for the aliases with default values.
        self.aliases = autonumerate(aliases)
              
    def __str__(self):
        """Represent as string.
        
        Overrides built-in list's representation.
        
        Returns:
            str
        """
        tmp = []
        for f, (myfile, myalias) in self.enum():
            tmp.append("\t".join([str(f), myfile, myalias]))
        tmp.append("")
        return "\n".join(tmp)
    
    def to_file(self, outfile=None, mode ='a'):
        """Save list as a text file that can be read back in.
        
        Args:
            outfile(str): Output file to write into. If omitted, it only 
                        returns the content as a print-ready string. 
                        (Default None)
            mode(str): Append ('a') or overwrite ('w'). (Default 'a')
        Returns:
            str: A print-ready multi-line string. This is returned even when an
                        output file is specified and written into.
        """
        result = ""
        for f, (myfile, myalias) in self.enum():
            result += myfile + "\t" + myalias + "\n"
        if outfile is not None:
            with open(outfile, mode) as out:
                out.write(result)
        return result    
        
    def enum(self):
        """Enumerate as (index, (filepath, filealias)).
        
        Returns:
            enumerator"""
        return enumerate(zip(self, self.aliases))
    
    def get(self, loc):
        """Access path and alias at specified location as tuple.
        
        Args:
            loc[int]: Index of item to get.
        Returns:
            (str,str): Tuple (path, alias).
        """
        return (self[loc], self.aliases[loc])
    
    def append(self, myfile, myalias=None):
        """Appends value to both the paths list and the aliases list.
        
        This method overrides the built-in append() of list. It is backwards
        compatible by automatically guessing an alias.
        This reduces the risk of the paths and aliases going out of sync due
        to items being manually added without updating the aliases.
        It is still possible to break the sync by manually adding items to the 
        aliases.
        
        Args:
            myfile(str): File (path will be expanded).
            myalias(str): Alias for the file (Default None).
        """
        if myfile is not None:
            myfile = os.path.abspath(os.path.expanduser(myfile))
        super(FilesList, self).append(myfile)
        if not myalias:
            myalias = self.autoalias(myfile)
        self.aliases.append(myalias)
        self.aliases = autonumerate(self.aliases)

    def populate_from_files(self, myfiles, colSep="\t"):
        """Parse the list of files from one or multiple text files. 
        
        Read in multiple lists of files from text and append them to the
        FilesList. All paths are automatically expanded and converted to 
        absolute paths. Because of this, each file may also be given a more 
        convenient alias. If no alias is given, the filename as supplied is 
        used as the alias instead. The paths are not checked for existence.
        
        Existing contents of the object are kept and the new contents are
        appended.
        
        Input file format (no spaces allowed inside names):
        
            #comment
            path1/file1     alias1-prefix    alias1-suffix1     alias1-suffix2
            path1/file2     alias1-prefix    alias1-suffix3
            path2/file3     alias3
            path3/file4
        
        Args:
            file[str]: A list of text files each containing a list of files.
            colSep(str): Column separator. (Default "\\t")
        Returns:
            FilesList: Returns self, to facilitate instantiation shortcuts.
        """
        # Read new list.
        paths = []
        for myfile in myfiles:
            with open(myfile, 'rU') as input:
                for line in input:
                    if line == "\n": 
                        # Skip empty lines.
                        continue
                    elif line[0] == '#':
                        # Skip comments.
                        continue
                    else:
                        fields = line.rstrip().split(colSep)
                        paths.append(fields[0])
                        # Store the alias for the file.
                        if len(fields) > 1:
                            self.aliases.append("_".join(fields[1:]))
                        # If an alias was not specified, re-use the filepath given.
                        else:
                            self.aliases.append(self.autoalias(fields[0]))
        # Expand to absolute paths and add to main self list.
        self.extend(expand_fpaths(paths))
        self.aliases = autonumerate(self.aliases)
        return self

    def populate_from_directories(self, dirpaths, patterns=None):
        """Get files based on naming patterns from within a list of directories.
        
        Useful for selecting among files that follow a naming convention. The
        convention is represented by a list of regex strings, at least one of
        which has to match.
        File paths will be expanded automatically. The filenames will be used 
        as the aliases.
        
        Existing contents of the object are kept and the new contents are
        appended.
        
        Args:
            dirpaths[str]: A list/FilesList of paths to directories from where
                        to get files.
            patterns[str]: A list of regex strings. Only files matching at least
                        one of these will be returned. The patterns will be 
                        matched anywhere in the filenames.
        Returns:
            FilesList: Returns self, to facilitate instantiation shortcuts.
        """
        rx = []
        if patterns:
            rx = [re.compile(p) for p in patterns]
        cwd = os.getcwd()
        for d in dirpaths:
            try:
                os.chdir(d)
                for f in os.listdir('./'):
                    if f in ["","\n",".",".."]: 
                        continue 
                    if not patterns:
                        # No filter.
                        self.append(f, self.autoalias(f))
                    else:
                        for p in rx:
                            if p.search(f):
                                self.append(f, self.autoalias(f))
                                break
            finally:
                # Ensure return to original directory to prevent errors in 
                # external code that uses the class.
                os.chdir(cwd)
        self.aliases = autonumerate(self.aliases)
        return self.sorted()
    
    # Helper function.
    @staticmethod
    def autoalias(pathname):
        """Strip a path to the base filename."""
        if pathname is None:
            return None
        else:
            return os.path.splitext(os.path.basename(pathname))[0]

    def sorted(self):
        """Sorted copy.
        
        Returns:
            FilesList
        """
        d = dict()
        for i, (myfile, myalias) in self.enum():
            d[myfile] = myalias 
        sk = natural_sorted(list(d.keys()))
        newFL = FilesList()
        for k in sk:
            newFL.append(k, d[k])
        return newFL

# Dispatch for similar functions.
_funcDispatch = {"dont_exist" :  dont_exist,
                 "arent_access" : arent_readable,
                 "are_empty" : are_empty,   # Doesn't make sense from command line 
                                            # because of lack of control of non-existence value.
                 "arent_text" : arent_text
                }


def store_metadata(flist, numoflines):
    """Store top lines of files into dictionary.
    
    Args:
        flist: A list or FilesList.
        numoflines(int): Number of lines to save.
    Returns:
        dict[]: The items of flist[] are used as keys.
    """
    metadata = dict()
    for myfile in flist:
        if myfile is None:
            fin = sys.stdin
        else:
            fin = open(myfile)
        lines = []
        for i in range(0, numoflines):
            lines.append(fin.readline())
        metadata[myfile] = "".join(lines)
        if fin != sys.stdin:
            fin.close()
    return metadata




#####   M A I N   #####


def main(args):
    """Provide command-line access to the module's functionality.

    The functionality and format of main is subject to change, as the module
    expands and evolves to suit my needs. Main() is not intended to be called 
    from within other code. 
    
    Optional short info is printed in commented lines. The info always 
    *succeeds* the relevant output, rather than precede it. This serves as 
    confirmation that the task completed. Calling details of the script are 
    recorded in commented lines at the top of the output.
    
    For more info on the functionality, read the above documentation of the 
    classes and functions. For usage syntax, execute the module with the -h 
    argument.
    
    
    """    
    # Organize arguments and usage help:
    parser = argparse.ArgumentParser(description="Provide INPUTTYPE and TARGETs \
     *before* providing any of the other parameters. This is due to many \
    parameters accepting an indefinite number of values. Only one task at a time.")
    
    # Input/Output.
    parser.add_argument('INPUTTYPE', type=str, choices=['L','T','D','P'],
                                help=" Specify the type of the TARGETs: \
                                'T' = The actual input filess. \
                                'L' = Text file(s) listing the input files. \
                                'P' = Get list of input files from STDIN pipe. \
                                'D' = Input data directly from STDIN pipe. \
                                ('D' is compatible with only some of the functions)")
    parser.add_argument('TARGET', type=str, nargs='*',
                                help=" The targets, space- or comma-separated. Usually files. \
                                Look into the specific task details below for special uses. \
                                Do not specify with INPUTTYPE 'P' or 'D'.")
    parser.add_argument('-O','--out', type=str, nargs=3,
                                help=" Send individual outputs to individual files instead of \
                                merging them to STDOUT. Output files will be like \
                                <out[0]>/<out[1]>target<out[2]>")
    # Parameters.
    parser.add_argument('-L','--log', action='store_true',
                                help=" Log this command to ./commands.log.")
    parser.add_argument('-c','--comments', action='store_true',
                                help=" Include commented info to STDOUT or files. (Default don't include)")
    parser.add_argument('-C','--STDERRcomments', action="store_false",
                                help=" Do NOT show info in STDERR. (Default show)")
    parser.add_argument('-s','--sep', type=str, default=["\t"], nargs='+',
                                help=" A list of input field separators. The first value \
                                will be used for all output. (Default \\t, bash syntax for tab: $'\\t').")
    parser.add_argument('-l','--labels', action='store_true',
                                help=" Discard column headers (first content line) in input files. (Default do not discard)")
    parser.add_argument('-r','--relabel', action='store_false',
                                help=" Do NOT create new column headers that reflect the origin of the columns. (Default create)")
    parser.add_argument('-i','--index', action='store_true',
                                help=" Use column 0 as row index. The index will always be included in the output. (Default no index)")
    parser.add_argument('-M','--metadata', type=int, default=0,
                                help=" Number of metadata lines at the \
                                beginning of input data (Default 0). Metadate will be read separately \
                                and re-added verbatim into the output.")
    # General tasks.
    parser.add_argument('--probe', type=str, choices=list(_funcDispatch.keys()),
                                help=" Do one of a these simple checks on the target files.")
    parser.add_argument('--dir', type=str, nargs='*',
                                help=" List the contents of the target paths. \
                                Full absolute file paths are returned. Each file is also given an alias. \
                                Supplying an optional list of regex patterns enables filtering of the result.")
    parser.add_argument('--link', type=str, nargs='+',
                                help=" Create symbolic links for the targets into the specified directory. \
                                Any additional values are used as respective names for the links, one for one, \
                                otherwise the aliases or basenames will be used, enumerated when necessary.")
    parser.add_argument('--loop', type=str, nargs='+',
                                help=" Repeat the specified shell command for each target value. \
                                The fist value of this parameter determines what the target values are: \
                                'S'= strings: paths/files/strings, 'R'= range: ranges of positive integers in x:y format. \
                                Target PLACEHOLDERS: ***full***, ***path***, ***file***, \
                                ***core***, ***alias***. \
                                If looping over a NUMERICAL RANGE use any of the last 3 placeholders. \
                                The nested command should be supplied as a list of components, not as a single string. \
                                Options intended for the nested command should be preceded \
                                by a '+' sign like this: '+-v'.")
    # Delimited file tasks.
    parser.add_argument('--swap', type=str,
                                help=" Replace all occurrences of the --sep values with the value supplied here.\
                                ** Bash syntax for tab: $'\\t'. Compatible with 'D' as INPUTTYPE.")
    parser.add_argument('--cntcols', action='store_true',
                                help="Count the number of fields in the first row of each target file.")
    parser.add_argument('--cols', nargs='+',
                                help="Extract the specified columns (named or 0-indexed) from each target. \
                                Column ranges in x:y format closed at both ends. \
                                Negative indices must be escaped first: \-1. Compatible with 'D' as INPUTTYPE.")
    parser.add_argument('--rndcols', type=int,
                                help="Randomly select this many columns from the target files. \
                                With --index, the index column will not be part of the random selection.")
    parser.add_argument('--appnd', action='store_true',
                                help="Append all the columns of the target files into a single table.")
    parser.add_argument('--valset', nargs=3,
                                help="Get the non-redundant set of values in the given row/column. \
                                Takes three arguments: (i) orientation 'r' for row or 'c' for column, \
                                (ii) position index of the row/column, (iii) repetition filter: \
                                'a' all values, 'u' unique values only, 'r' only values with two or more instances.")
    params = parser.parse_args(args)
    
    # INPUT ###################################################################
    
    targets = []
    for t in params.TARGET:
        v = t.split(",")
        if len(v) == 1:
            targets.append(t)
        else:
            targets.extend(v)
    flist = None
    if params.INPUTTYPE == 'P':
        # Read files list from STDIN
        flist = FilesList()
        for line in sys.stdin:
            fields = line.rstrip("\n").split("\t")
            if fields[0] != "":
                try:
                    flist.append(fields[0], fields[1])
                except IndexError:
                    flist.append(fields[0])
    elif params.INPUTTYPE == 'L':
        # Create the FilesList, by appending the contents of all provided lists.
        flist = FilesList().populate_from_files(targets)
    elif params.INPUTTYPE == 'T':
        # Create the FilesList by supplying a direct list of files.
        flist = FilesList(targets)
    elif params.INPUTTYPE == 'D':
        # Data will be read from STDIN. No files needed. Make an empty list. 
        # Not all functions will switch to STDIN given this. Several will simply do nothing.
        flist = FilesList()
    else:
        sys.exit(ml.errstring("Unknown INPUTTYPE."))
        
    # Metadata. ---------------------------------------------------------------
    metadata = ""
    if params.metadata:
        metadata = store_metadata(flist, params.metadata)
    
    # OUTPUT ##################################################################
    
    outdir, outpref, outsuff = None, None, None
    if params.out:
        outdir = expand_fpaths([params.out[0]])[0]
        outpref = params.out[1]
        outsuff = params.out[2]
    
    # CALL DETAILS ############################################################
    
    if params.log:
        ml.log_command()
    if params.STDERRcomments:
        sys.stderr.write(ml.paramstring())
    
    # TASKS ###################################################################
     
    # Simple file PROBEs. -----------------------------------------------------
    if params.probe:
        result = _funcDispatch[params.probe](flist)
        try:
            if params.comments:
                sys.stdout.write(ml.paramstring())
            print(result.to_file())
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("probing"))
        except IOError:
            pass
    
    
    # Filter DIRECTORY contents. ----------------------------------------------
    elif params.dir is not None:
        result = FilesList().populate_from_directories(flist, params.dir)
        try:
            if params.comments:
                sys.stdout.write(ml.paramstring())
            sys.stdout.write(result.to_file())
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("listing"))
        except IOError:
            pass
         
         
    # LOOP arbitrary command. -------------------------------------------------
    elif params.loop:
        if params.loop[0] == 'R':
            # Generate the range.
            myrange = []
            for t in targets:   # Look for multiple ranges.
                v = t.split(":")
                if len(v) > 1:
                    myrange.extend(list(range(int(v[0]), int(v[1]) + 1)))
                else:
                    sys.exit(ml.errstring("No numeric ranges specified. Use -h for help with the newest syntax."))
            flist = FilesList(myrange)
        # Strip left and/or right padding first.
        command = []
        for c in params.loop[1:]:
            command.append(c.lstrip("+"))
        try:
            do_foreach(flist, command, out=(outdir, outpref, outsuff),
                       progress=(params.STDERRcomments), comments=params.comments,
                       log=params.log)
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("looping-"+ params.loop[0]))
        except IOError:
            pass
    
    
    # Symbolic LINKS. ---------------------------------------------------------
    elif params.link:
        slink(flist, dir=params.link[0], aliases=params.link[1:])
        if params.STDERRcomments:
            sys.stderr.write(ml.donestring("linking"))
             
             
    # SWAP substrings. --------------------------------------------------------
    elif params.swap is not None:
        result = swap_strFiles(flist, insep=params.sep, outsep=params.swap)
        # Create output filenames, if applicable. If [], then STDOUT.
        outfiles = make_names(flist.aliases, (outdir, outpref, outsuff))
        outstream = sys.stdout
        # I need the for loop to iterate at least once. Relevant for STDIN input, since I have no input files listed then.
        if flist == []:
            flist.append("<STDIN>")
        # Print the converted data.
        for i, (myfile, myalias) in flist.enum():
            if outfiles:
                # Send to individual file instead of STDOUT.
                outstream = open(outfiles[i], 'w')
            try:
                if params.comments:
                    # Embed call info at beginning of output. More useful there when outputting to files.
                    outstream.write(ml.paramstring("SOURCE: " + myfile))
                outstream.write(result[i].rstrip("\n") +"\n")
            except IOError:
                pass
            finally:
                if outfiles:
                    # Don't want to accidentally close STDOUT.
                    outstream.close()
        if params.STDERRcomments:
            try:
                sys.stderr.write(ml.donestring("swapping delimiters"))
            except IOError:
                pass
    
         
    # Get COLUMNS or RANDOM columns. (most code shared) -----------------------
    elif params.cols or params.rndcols:
        # Create output filenames, if applicable. If [], then STDOUT.
        outfiles = make_names(flist.aliases, (outdir, outpref, outsuff))
        outstream = sys.stdout
        merge = False if outfiles else True
        # Determine if using index, and assign appropriate value.
        idx = None
        if params.index:
            idx = 0
        else:
            idx = None
        # Extract data.
        result = None
        if params.cols:
            cols = []
            for p in params.cols:   # space separated arguments
                cols.extend(p.split(","))  # comma separated arguments
            # Get the specified columns.
            result = get_columns(flist, cols=cols, colSep=params.sep, 
                                header=params.labels, merge=merge, index=idx)
        else:
            # Get random columns.
            result = get_random_columns(flist, k=params.rndcols, colSep=params.sep, 
                                        header=params.labels, merge=merge, index=idx)
        # I need the for loop to iterate at least once. Relevant for STDIN input, since I have no input files listed then.
        if flist == []:
            flist.append("<STDIN>")
        if merge:
            try:
                if params.comments:
                    # Embed call info at beginning of output.
                    outstream.write(ml.paramstring("SOURCE: " + myfile))
                if params.metadata:
                    # Dump all the metadata from all the merged input sources.
                    for i, (myfile, myalias) in flist.enum():
                        outstream.write(metadata[myfile])
                outstream.write( result[0].to_csv(header=params.relabel, index=params.index, sep=params.sep[0]))
            except IOError:
                pass
        else:    
            for i, (myfile, myalias) in flist.enum():
                outstream = open(outfiles[i], 'w')
                try:
                    if params.comments:
                        # Embed call info at beginning of output.
                        outstream.write(ml.paramstring("SOURCE: " + myfile))
                    if params.metadata:
                        outstream.write(metadata[myfile])
                    outstream.write( result[i].to_csv(header=params.relabel, index=params.index, sep=params.sep[0]))
                except IOError:
                    pass
                finally:
                    outstream.close()
        if params.STDERRcomments:
            try:
                if params.cols:
                    sys.stderr.write(ml.donestring("getting columns, index "+ str(idx is not None)))
                else:
                    sys.stderr.write(ml.donestring("getting random columns, index "+ str(idx is not None)))
            except IOError:
                pass
        
    
    # APPEND columns. ---------------------------------------------------------
    elif params.appnd:
        idx = None
        if params.index:
            idx = 0
        df = append_columns(flist, colSep=params.sep, header=params.labels, index=idx)
        try:
            if params.comments:
                ml.parastring()
            if params.metadata:
                # Dump all the metadata from all the merged input sources.
                for i, (myfile, myalias) in flist.enum():
                    outstream.write(metadata[myfile])
            sys.stdout.write(df.to_csv(sep=params.sep[0], header=params.relabel, index=params.index))
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("appending columns, index "+ str(idx is not None)))
        except IOError:
            pass
 
 
    # COUNT columns. ----------------------------------------------------------
    elif params.cntcols:
        result = count_columns(flist, params.sep)
        try:
            if params.comments:
                sys.stdout.write(ml.paramstring())
            for f, (myfile, myalias) in flist.enum():
                print("\t".join([str(result[f]), myalias, myfile]))
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("counting columns"))
        except IOError:
            pass
        
    
    # SET of values in row/column. --------------------------------------------
    elif params.valset:
        nest = get_valuesSet(flist, axis=params.valset[0], index=params.valset[1], filter=params.valset[2], colSep=params.sep)
        try:
            if params.comments:
                sys.stdout.write(ml.paramstring())
            for f, (myfile, myalias) in flist.enum():
                print("".join([myfile, "\t", str(nest[f])]))
            if params.STDERRcomments:
                sys.stderr.write(ml.donestring("obtaining set of values."))
        except IOError:
            pass
        
        
     
#     # All done.     
#     if params.STDERRcomments:
#         sys.stderr.write(ml.donestring())   
        
        
        


#####    E X E C U T I O N   #####


# Call main only if the module was executed directly.
if __name__ == "__main__":
    main(sys.argv[1:])
    
    sys.exit(0)
    
#EOF