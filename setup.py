import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# Required package

REQUIRED = ['pandas', 'numpy', 'py4cytoscape',
            'seaborn', 'more-itertools', 'argparse']

# This call to setup() does all the work
setup(
    name="metabodirect",
    version="0.3.2",
    description="Analyze FT-ICR-MS data",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/Coayala/MetaboDirect",
    author="Christian Ayala",
    author_email="cayalaortiz@email.arizona.edu",
    license="LICENSE",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    packages=["metabodirect"],
    include_package_data=True,
    install_requires=REQUIRED,
    entry_points={
        "console_scripts": [
            "metabodirect=metabodirect.__main__:main",
            "test_normalization=metabodirect.test_normalization:main",
            "create_networks=metabodirect.create_networks:main"
        ]
    },
)
