from setuptools import setup, find_packages
#from Cython.Build import cythonize
from distutils.extension import Extension

with open("README.md","r") as fh:
	long_description = fh.read()

try:
    from Cython.Distutils import build_ext
except ImportError:
	ext_modules = [
        Extension("igblast_parser.module1", ["compiled/module1.c"]),
    ]
	cmdclass = {}
else:
	#ext_modules = cythonize("compiled/module1.pyx")
	cmdclass = {'build_ext': build_ext}
	ext_modules = [
        Extension("igblast_parser.module1", ["compiled/module1.pyx"]),
    ]

setup(
	name = "igblast_parser",
	version = "0.0.3",
	author = "Ariel Erijman & Brd Langhorst",
	author_email = "aerijman@fredhutch.org",
	description = "Parser of Igblast results into a csv file",
	long_description = long_description,
	long_description_content_type = "text/markdown",
	url = "https://github.com/aerijman/igblast-parser",
	packages = find_packages(),
	#classifiers = [
	#	"Programming language ::  Python :: 3",
	#	"Licence :: OSI approved :: MIT Licence",
	#	"Operating system :: OS Independent",
	#],
	python_requires='>=3.6',
	install_requires = ['pandas', 'numpy'],
	
	ext_modules = ext_modules,
	cmdclass = cmdclass,
	scripts = ['bin/igblast-parser'],
)
