from setuptools import setup, find_packages

setup(name='ropeat',
        version='v0.0.0',
        packages=find_packages(),
        package_data={'ropeat/models': ['*.pkl']},
        include_package_data=True)
