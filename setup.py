import sys

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

setup(
    name="flows1d0d3d",
    version="0.0.1",
    description="Package for calculating the flow and transport inside a vascular network",
    author="Andreas Wagner, Prashant K. Jha, Tobias Koeppl, Marvin Fritz",
    license="boost software license 1.0",
    packages=find_packages(where="python/src"),
    package_dir={"": "python/src"},
    cmake_install_dir="python/src/flows1d0d3d",
    cmake_args=['-DLibMacrocirculation_Enable_Apps:BOOL=OFF',
                '-DLibMacrocirculation_Enable_Tests:Bool=OFF',
                '-DLibMacrocirculation_Enable_Python_Bindings=ON'],
    include_package_data=True,
    extras_require={"test": ["pytest"]},
    python_requires=">=3.6",
)