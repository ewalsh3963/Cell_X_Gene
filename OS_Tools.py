import os
import glob
import shutil

def ensure_directory(path, critical=False):
    if not os.path.isdir(path):
        if critical:
            raise FileNotFoundError(f"Critical directory does not exist: {path}")
        else:
            os.makedirs(path)

def find_files(parent_dir, extension, wd=None, recursive=False):
    """
    Searches for files within a specified directory using glob patterns.

    Args:
        parent_directory (str): The starting directory for the search.
                                This can be an absolute or relative path.
        extension (str): The file extension to search for (e.g., '.txt', '.py', 'jpg').
                         Include the dot if you want to match files with an extension.
                         Use an empty string '' or '*' to match any extension (if combined with *).
        working_directory (str, optional): An optional working directory to change into
                                           before performing the search. If None, the
                                           current working directory is used.
        recursive (bool): If True, the search will include subdirectories recursively.
                          If False, only the parent_directory is searched.

    Returns:
        list: A list of paths to the found files, relative to the working_directory
              or the directory where the script is run if no working_directory is specified.
    """
    original_cwd = os.getcwd()  # Store the original current working directory
    found_files = []

    try:
        if wd:
            os.chdir(wd)
            print(f"Changed working directory to: {os.getcwd()}")
        else:
            print(f"Searching from current working directory: {os.getcwd()}")

        # Ensure the extension starts with a dot if it's not empty and doesn't already have one
        if extension and not extension.startswith('.'):
            search_extension = f"*.{extension}"
        elif extension == "": # If extension is empty, match any file name (e.g., "file." for files with no extension)
            search_extension = "*" # Match any file
        else: # e.g. '.txt' or '*'
            search_extension = f"*{extension}"

        # Construct the glob pattern
        if recursive:
            # Use '**' to match any directories and subdirectories recursively
            # os.path.join handles path separators correctly across OS
            search_pattern = os.path.join(parent_dir, '**', search_extension)
            found_files = glob.glob(search_pattern, recursive=True)
        else:
            # Search only in the parent_directory
            search_pattern = os.path.join(parent_dir, search_extension)
            found_files = glob.glob(search_pattern, recursive=False)

    except FileNotFoundError:
        print(f"Error: Directory '{wd or os.getcwd()}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        # Always change back to the original working directory
        if os.getcwd() != original_cwd:
            os.chdir(original_cwd)
            print(f"Changed back to original working directory: {os.getcwd()}")

    return found_files

def check_file(file, critical=False, rewrite=True):
    """
    Checks if a file exists and initializes it if rewrite is True or it doesn't exist.

    Parameters:
    - file (str): Path to the file.
    - critical (bool): If True, raises an error if the file does not exist.
    - rewrite (bool): If True and file exists, overwrite allowed (file will be truncated).
                      If False and file exists, do not rewrite.

    Returns:
    - bool: True if file is ready (exists or created),
            False if file exists and rewrite is False.
    """
    file_exists = os.path.exists(file)
    
    if file_exists:
        if rewrite:
            # Initialize (truncate) the existing file
            with open(file, 'w'):
                pass
            return True
        else:
            # File exists and rewrite not allowed
            return False
    else:
        if critical:
            raise FileNotFoundError(f"Critical: file '{file}' does not exist.")
        # Create (initialize) the file
        with open(file, 'w'):
            pass
        return True

## Writing log files 
class Logs:
    def __init__(self, log_file):
        self.log_file = log_file

    def append_logs(self, text):
        """Append the given text to the log file, adding a newline if needed."""
        with open(self.log_file, 'a', encoding='utf-8') as f:
            if not text.endswith('\n'):
                text += '\n'
            f.write(text)

def clear_directory(path):
    """
    Recursively deletes all contents inside the given directory,
    but leaves the directory itself intact.

    Parameters:
    - path (str): Path to the target directory.
    """
    if not os.path.isdir(path):
        return
    
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        try:
            if os.path.isfile(item_path) or os.path.islink(item_path):
                os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
        except Exception as e:
            print(f"Failed to delete {item_path}: {e}")

def delete_file_if_exists(filepath):
    """
    Deletes the file at the specified path if it exists.

    Parameters:
    - filepath (str): The path to the file.

    Returns:
    - bool: True if file existed and was deleted, False if file did not exist.
    """
    if os.path.isfile(filepath):
        os.remove(filepath)
        print(f"Deleted file: {filepath}")
        return True
    else:
        print(f"File not found: {filepath}")
        return False