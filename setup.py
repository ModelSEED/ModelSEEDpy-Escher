from setuptools import setup, find_packages

setup(name='modelseedpy_escher',
      version='0.2.0',
      description='Escher viewer for ModelSEED',
      url='https://github.com/ModelSEED/modelseed-escher',
      author='Filipe Liu',
      author_email='fliu@anl.gov',
      license='MIT',
      
      packages=find_packages(),
      install_requires=[
          "cobra >= 0.15.3",
          "cobrakbase >= 0.1.9",
          "escher == 1.6.0",
      ],
      zip_safe=False)