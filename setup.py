from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name = "pyphot",
    version = 1.2.1,
    description = "A tool for computing photometry from spectra",
    long_description = readme(),
    author = "Morgan Fouesneau",
    author_email = "",
    url = "https://github.com/mfouesneau/pyphot",
    packages = find_packages(),
    package_data = {'pyphot':['libs/*'], 
                    'pyphot.ezunits':['default_en.txt']},
    include_package_data = True,
    classifiers=[
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=["numpy", "scipy", "matplotlib","tables", "pandas", "astropy" ]
)
