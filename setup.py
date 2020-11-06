from setuptools import find_packages, setup

from pathlib import Path

with open("README.md", "r") as fh:
    long_description = fh.read()

def parse_requirements(requirements_path):
    with open(Path(__file__).parent / requirements_path) as f:
        return f.read().splitlines()

requirements = parse_requirements("requirements.txt")
print("REQUIREMENTS: ", requirements)

setup(
    name="scEasyMode",
    version="1.0.1",
    author="Johnny Yu",
    author_email="johnny.yu@ucsf.edu",
    description="Wrappers for automating single cell workflows in python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/johnnyUCSF/scEasyMode",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
