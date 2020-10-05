import setuptools

with open("README.md", "r") as readme:
    long_description = readme.read()

with open("requirements.txt", "r") as f:
    install_requires = f.read()

setuptools.setup(
    name="virtool_cli",
    version="0.0.1",
    description="a command line interface for working with Virtool data",
    long_description=long_description,
    license="MIT",
    platforms="linux",
    python_requires=">=3.6",
    url="https://github.com/virtool/virtool-cli",
    packages=setuptools.find_packages(exclude="tests"),
    install_requires=install_requires,
    entry_points={
        "console_scripts": [
            "virtool = virtool_cli.run:cli",
        ],
    }
)
