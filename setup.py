from distutils.core import setup

setup(
    name='Contrast',
    author="Enrique Paillas",
    author_email="enrique.paillas@uwaterloo.ca",
    version='0.1dev',
    url='https://github.com/epaillas/contrast',
    packages=['contrast'],
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),

    install_requires=[
                    'julia'
                    'numpy',
                    'jupyter',
                    'matplotlib',
                    ],

    zip_safe=False
)
