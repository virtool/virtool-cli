import setuptools


setuptools.setup(
    name="virtool_cli",
    packages=setuptools.find_packages(exclude="tests"),
    version="0.0.1",
    description="a command line interface for working with Virtool data",
    install_requires=["click"],
    url="https://github.com/virtool/virtool-cli",
    entry_points={
        "console_scripts": [
            "virtool = virtool_cli.run:cli",
        ],
    }
)
