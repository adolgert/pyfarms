from setuptools import setup

setup(name="pyfarms",
    version="0.1",
    description="Farm disease spread model",
    long_description="""An equivalent to NAADSM.
    """,
    classifiers=[
      "Development Status :: 3 - Alpha",
      "Intended Audience :: Science/Research",
      "Natural Language :: English",
      "License :: OSI Approved :: BSD License",
      "Programming Language :: Python :: 3.4",
      "Topic :: Scientific/Engineering :: Mathematics"
    ],
    keywords="stochastic dynamic simulation Markov Gillespie FMD HPAI",
    url="http://github.com/adolgert/pyfarms",
    author="Drew Dolgert",
    author_email="ajd27@cornell.edu",
    license="BSD 3-clause",
    packages=["pyfarms"],
    install_requires=[
        "numpy",
        "scipy",
        "gspn"
    ],
    include_package_date=True,
    zip_safe=False,
    entry_points={
      "console_scripts" : ["farmspread=pyfarms.command_line:main"]
    })
