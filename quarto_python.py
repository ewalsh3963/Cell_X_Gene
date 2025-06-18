import subprocess
import os
import functools
import shutil
import pdb

def ensure_same_dir_after_execution(f):

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        cwd = os.getcwd()
        try:
            os.chdir(os.path.dirname(__file__))
            return f(*args, **kwargs)
        finally:
            os.chdir(cwd)

    return wrapper


@ensure_same_dir_after_execution
def call_render_shell_script(set_path=False, ):
    script_path = 'render'
    if script_path not in os.listdir():
        os.chdir(os.path.dirname(__file__))
    if script_path not in os.listdir():
        raise FileNotFoundError(f"Could not find {script_path}")
    process = subprocess.run(
        [script_path],
        shell=True,
        capture_output=True,
        # stdout=subprocess.PIPE,
        # stderr=subprocess.PIPE,
    )
    print(process.stdout.decode())
    print(process.stderr.decode())


def get_quarto_path():
    path_env = os.getenv("QUARTO_PATH")
    if path_env is None:
        return shutil.which("quarto")
    else:
        return path_env


# can pass params as dict e.g. { 'data_path' : 'data/study_1.csv'}
def quarto_direct(
    *args,
    quarto_cmd='render',  # can be 'preview'
    quarto_execute_params_as_dict=None,
    clean_site_libs_after_render=True,
    **kwargs,
):
    quarto_executable = get_quarto_path()
    cmd = [
        quarto_executable,
        quarto_cmd,
        *flatten((f'--{key}', f'{value}') for key, value in kwargs.items()),
        *args,
    ]
    if quarto_execute_params_as_dict is not None:
        cmd.extend(
            [f'-P {k}:{v}' for k, v in quarto_execute_params_as_dict.items()])
    print(' '.join(cmd))
    
    try:
        process = subprocess.Popen(cmd)
        process.wait()
    except KeyboardInterrupt:
        # this makes it possible to use the 'preview' command and stop it gracefully
        print("KILLING PROCESS")
        process.kill()

    if clean_site_libs_after_render:
        clean_site_libs()


def flatten(l):
    return [item for sublist in l for item in sublist]


def clean_site_libs(dir='.', ):
    cmd = f"find {dir} -type d -name 'site_libs' -exec rm -rf {{}} +"
    process = subprocess.run(cmd, shell=True)
    print("Cleaned up site_libs")


# TODO: possibly add an improved interface over Quarto's for getting parameters associated with a profile
# (i.e. parse YAML here)
def quarto_profile(
    profile_name,
    quarto_cmd='render',  # can be 'preview'
    use_profile_config_for_parameters=True,
    clean_site_libs_after_render=True,
    **kwargs,
):
    cmd = [
        f'--profile',
        profile_name,
    ]
    if use_profile_config_for_parameters:
        cmd.extend([
            f'--execute-params',
            f'_quarto-{profile_name}.yml',
        ])
    return quarto_direct(
        *cmd,
        quarto_cmd=quarto_cmd,
        clean_site_libs_after_render=clean_site_libs_after_render,
        **kwargs,
    )


# This does what the shell `render` script does, but from Python
# TODO: possibly pass parameters or add YAML to populate navbar for each profile and its `output-dir`` (specified in _quarto-<profile>.yml)
def render_all_profiles():
    profiles = [
        file.removeprefix('_quarto-').removesuffix('.yml')
        for file in os.listdir()
        if file.startswith('_quarto-') and file.endswith('.yml')
    ]
    for profile in profiles:
        quarto_profile(
            profile,
            use_profile_config_for_parameters=True,
            clean_site_libs_after_render=True,
        )


def preview_profile(profile_name, **kwargs,):
    quarto_profile(
        profile_name,
        quarto_cmd='preview',
        # port=8080,  # does not work for some reason when called from subprocess
        **kwargs,
    )


def test_call_render_shell_script():
    call_render_shell_script()


if __name__ == "__main__":
    render_all_profiles()
    # quarto_profile('study_1')
    # preview_profile(profile_name='study_1')
