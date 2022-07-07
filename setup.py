import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='tierra',
     version='1.0',
     author="Prajwal Niraula",
     author_email="pniraula@mit.edu",

     description="1D transmission spectroscopy code.",
     long_description=long_description,

     long_description_content_type="text/markdown",
     url="https://github.com/disruptiveplanets/tierra",
     packages=setuptools.find_packages(),

     classifiers=[
         "Programming Language :: Python :: >3.5",
         "License :: MIT License",
         "Operating System :: LINUX",
     ],

 )
