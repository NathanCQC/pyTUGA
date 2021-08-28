import os
import sys

from setuptools.command.bdist_egg import bdist_egg as _bdist_egg
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.install import install as _install
from setuptools.command.build_ext import build_ext

from distutils.command.build import build as _build
from distutils.command.bdist import bdist as _bdist
from distutils.command.sdist import sdist as _sdist

from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
from Cython.Distutils import Extension

CUSTOM_OPTION = ("pure", None, "Build pure python distribution")

def _check_platform():
    if os.name in ("nt", "dos"):
        msg = "pyTUGA does not support Windows"
        sys.exit(msg)


def _check_python():
    v = sys.version_info
    if v[:2] < (3, 7):
        error = "ERROR: pyTUGA requires Python version 3.7, 3.8 or 3.9."
        sys.exit(error)


def get_version():
    version = {}
    with open("pyTUGA/version.py") as fp:
        exec(fp.read(), version)
    return version


class build(_build):
    user_options = _build.user_options + [CUSTOM_OPTION]
    boolean_option = _build.boolean_options + [CUSTOM_OPTION[0]]

    def initialize_options(self):
        _build.initialize_options(self)
        self.pure = None

    def finalize_options(self):
        _build.finalize_options(self)


class install(_install):
    user_options = _install.user_options + [CUSTOM_OPTION]
    boolean_options = _install.boolean_options + [CUSTOM_OPTION[0]]

    def initialize_options(self):
        _install.initialize_options(self)
        self.pure = None

    def finalize_options(self):
        _install.finalize_options(self)
        build_cmd = self.get_finalized_command("build")
        if build_cmd.pure is None:
            build_cmd.pure = self.pure

    def run(self):
        """Running grandparent install.run()"""
        return super(_install, self).run()


class bdist(_bdist):
    user_options = _bdist.user_options + [CUSTOM_OPTION]
    boolean_option = _bdist.boolean_options + [CUSTOM_OPTION[0]]

    def initialize_options(self):
        _bdist.initialize_options(self)
        self.pure = None

    def finalize_options(self):
        _bdist.finalize_options(self)
        build_cmd = self.get_finalized_command("build")
        if build_cmd.pure is None:
            build_cmd.pure = self.pure


class bdist_wheel(_bdist_wheel):
    user_options = _bdist_wheel.user_options + [CUSTOM_OPTION]
    boolean_options = _bdist_wheel.boolean_options + [CUSTOM_OPTION[0]]

    def initialize_options(self):
        _bdist_wheel.initialize_options(self)
        self.pure = None

    def finalize_options(self):
        build_cmd = self.get_finalized_command("build")
        if build_cmd.pure is None:
            build_cmd.pure = self.pure
        _bdist_wheel.finalize_options(self)
        if not self.pure:
            self.compile = self.optimize = False
            self.root_is_pure = False


class build_py(_build_py):
    user_options = _build_py.user_options + [CUSTOM_OPTION]
    boolean_option = _build_py.boolean_options + [CUSTOM_OPTION[0]]

    def initialize_options(self):
        _build_py.initialize_options(self)
        self.pure = None

    def finalize_options(self):
        _build_py.finalize_options(self)
        self.set_undefined_options("build", ("pure", "pure"))

    def is_compile(self, module_name):
        """Don't compile __init__.py files."""
        if self.pure:
            return False
        else:
            return module_name.find("__init__") == -1

    def run(self):
        if not self.pure:
            if not self.distribution.ext_modules:
                self.distribution.ext_modules = []
            for (pkg, mod, pth) in self.find_all_modules():
                if self.is_compile(mod) and mod != "version":
                    self.distribution.ext_modules.append(
                        Extension(
                            ".".join([pkg, mod]),
                            [pth],
                            cython_c_in_temp=True,
                            cython_directives={"language_level": 3},
                        )
                    )
        return _build_py.run(self)

    def build_packages(self):
        for package in self.packages:
            package_dir = self.get_package_dir(package)
            modules = self.find_package_modules(package, package_dir)
            for (package_, module, module_file) in modules:
                assert package == package_
                if not self.is_compile(module) or module == "version":
                    _build_py.build_module(self, module, module_file, package)


class bdist_egg(_bdist_egg):
    """
    Hacked bdist_egg to always build a pure python egg.
    Needed to generate a *.pyc-only byte-compiled python
    wheel without *.py source.
    """

    def initialize_options(self):
        _bdist_egg.initialize_options(self)
        self.pure = None

    def finalize_options(self):
        _bdist_egg.finalize_options(self)
        build_py_cmd = self.get_finalized_command("build_py")
        if build_py_cmd.pure is None:
            build_py_cmd.pure = True


class disabled_bdist_egg(_bdist_egg):
    """
    Disabled version of bdist_egg. Prevents 'python setup.py install' running
    setuptools' default easy_install command. We do not want to run easy_install.
    """

    def run(self):
        msg = r"""
        Aborting implicit building of eggs. 
        Use `pip install .` to install from source.
        """
        sys.exit(msg)


class sdist(_sdist):
    """
    Disabled Version of sdist. Prevents 'python setup.py sdist' from creating a pure python
    source distribution and packaging as a tarbell, zip file etc
    """

    def run(self):
        msg = r"""
        Aborting source distribution build.
        Use 'python setup.py bdist' to
        create a cythonized binary distribution.
        """
        sys.exit(msg)


if __name__ == "__main__":
    from setuptools import setup

    _configuration = {
        "cmdclass": {
            "install": install,
            "bdist": bdist,
            "bdist_wheel": bdist_wheel,
            "build": build,
            "build_py": build_py,
            "build_ext": build_ext,
            "bdist_egg": bdist_egg if "bdist_egg" in sys.argv else disabled_bdist_egg,
        },
        "use_scm_version": {
            "root": ".",
            "relative_to": __file__,
            "version_scheme": "release-branch-semver",
            "write_to": "pyTUGA/version.py",
            "git_describe_command": 'git describe --dirty --tags --long --match "*[0-9]*" --abbrev=14',
        },
        "setup_requires": ["setuptools_scm"],
    }

    _check_platform()
    _check_python()
    setup(**_configuration)
