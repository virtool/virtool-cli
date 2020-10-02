import setuptools


setuptools.setup(
    name="virtool-cli",
    packages=setuptools.find_packages(),
    version="0.0.1",
    description="a command line interface for working with Virtool data",
    py_modules=["run"],
    install_requires=['Click'],
    url="https://github.com/virtool/virtool-cli",
    entry_points={
        'console.scripts': [
            'virtool = run:cli',
        ],
    }
)
