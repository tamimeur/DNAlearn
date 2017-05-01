from setuptools import setup

setup(name='optseq',
      version='0.1',
      description='Recommender system for sequence design and optimization',
      url='https://github.com/tamimeur/OptSeq.git',
      author='Leli Ami',
      author_email='tamimeur@gmail.com',
      license='MIT',
      packages=['optseq'],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'optseq = optseq.cli:main',
        ],
      })
