from setuptools import setup, find_packages

def get_readme():

    with open('README.md') as f:
        return f.read()

setup(

  name='sobols',
  description='Variance-based sensitivity analysis.',
  long_description=get_readme(),
  version='0.0.0',
  license='MIT License',
  author='Antonio Mastromarino, Sara Sottile, Mattia Sensi',
  author_email='antonio.mastromarino@wolfson.ox.ac.uk',
  maintainer='Antonio Mastromarino, Sara Sottile, Mattia Sensi',
  maintainer_email='antonio.mastromarino@wolfson.ox.ac.uk',
  url='https://github.com/AntonioVitoMastromarino/SobolS',

  packages=find_packages(include=('sobols', 'sobols.*')),
    install_requires=[
      'numpy',
    ],
    extras_require={
      'docs': [
        'sphinx',
        'sphinx_rtd_theme',
      ],
    },
)
