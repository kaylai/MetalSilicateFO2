import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fO2calculate",
    version="0.1.0",
    author="Kayla Iacovino",
    author_email="kaylaiacovino@gmail.com",
    description=("A python library for calculating the oxygen fugacity of metal-silicate "
                 "systems using non-ideal activities for Fe-FeO and Si-SiO2 equilibria."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kaylai/MetalSilicateFO2",
    packages=setuptools.find_packages(),
    install_requires=[
            'pandas',
            'numpy',
            'mendeleev'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
