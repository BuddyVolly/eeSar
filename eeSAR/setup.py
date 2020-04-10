import os
from setuptools import setup, find_packages


def parse_requirements(file):
    return sorted(set(
        line.partition('#')[0].strip()
        for line in open(os.path.join(os.path.dirname(__file__), file))
    )
                  -set('')
                  )


setup(
    name='eeSar',
    packages=find_packages(),
    include_package_data=True,
    version='0.0.1',
    description='Helper functions for SAR processsing with GEE',
    install_requires=parse_requirements('requirements.txt'),
    url='https://github.com/BuddyVolly/eeSAR',
    author='Andreas Vollrath',
    author_email='andreas.vollrath[at]esa.int',
    license='MIT License',
    keywords=['Sentinel-1', 'ESA', 'SAR', 'Radar',
              'Earth Observation', 'Remote Sensing',
              'Synthetic Aperture Radar', 'Google Earth Engine'],
    zip_safe=False,
)