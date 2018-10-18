import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
	name="targetDB",
	version="1.1.4",
	author="Stephane De Cesco",
	author_email="stephane.decesco@ndm.ox.ac.uk",
	description="Package with an application to generate report on potential drug targets",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/sdecesco/targetDB",
	packages=setuptools.find_packages(),
	package_data={'targetDB': ['data/*.zip','ml_data/*.zip', 'LICENSE']}
	, classifiers=("Programming Language :: Python :: 3", "Development Status :: 4 - Beta", "Intended Audience :: Science/Research","License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
		"Operating System :: OS Independent", "Programming Language :: Python :: 3 :: Only"
	),
	install_requires=['biopython', 'scipy', 'matplotlib', 'pandas', 'intermine', 'opentargets', 'xmltodict'],
	entry_points={'console_scripts': ['target_DB=targetDB.druggability_DB:entry_point',
	                                  'target_REPORT=targetDB.druggability_report:entry_point']}
)
