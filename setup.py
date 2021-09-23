from setuptools import setup, find_packages

setup(name='ModelSEEDpy-Escher',
      version='0.2.2',
      description='Escher viewer for ModelSEED',
      url='https://github.com/ModelSEED/ModelSEEDpy-Escher',
      author='Filipe Liu',
      author_email='fliu@anl.gov',
      license='MIT',
      
      packages=find_packages(),
      install_requires=[
          "networkx >= 2.4",
          "cobra >= 0.15.3",
          "cobrakbase >= 0.1.9",
          "escher == 1.6.0",
      ],
      zip_safe=False)
